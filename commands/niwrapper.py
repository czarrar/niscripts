#!/usr/bin/env python

import argparse, os, yaml, sys
import os.path as op
from datetime import datetime
from analysis.base import SubjectBase
from subprocess import Popen

config_file = "niwrapper.yaml"

def die(msg):
    print "ERROR: %s" % msg
    SystemExit(2)

class NiWrapper(SubjectBase):
    def __init__(self, config_file):
        # Load config
        if not op.isfile(config_file):
            raise Exception("Cannot find configuration file: %s" % config_file)
        f = file(config_file, 'r')
        self.config = yaml.load(f)
        f.close()
        
        # Get context
        template_context = self.config.pop("vars", {})
        
        # Setup
        self._parser_help = {}
        self._commands_opts = {}
        for k,opts in self.config.iteritems():
            # help
            if k in self._parser_help:
                die("Duplicate sub-command %s" % k)
            self._parser_help[k] = opts.pop("help", "sub-command")
            # command
            if k in self._commands_opts:
                die("Duplicate sub-command %s" % k)
            try:
                p = opts.pop("program")
            except KeyError:
                die("Must specify program in config file")
            self._commands_opts[k] = (p, opts)
        
        self._is_parsed = False
        return
    
    def setup(self, run_keys, log_dir, subjects, sge, sge_opts, verbosity, dry_run=False, 
                     **cmds):
        # Set subjects
        new_subjects = []
        for s in subjects:
            if s[0] == "@":
                f = file(s, 'r')
                tmps = f.readlines()
                tmps = [ l.strip() for l in tmps if l.strip() ]
                new_subjects.extend(tmps)
            new_subjects.append(s)
        subjects = new_subjects
        
        # Save
        self.run_keys = run_keys
        self.subjects = subjects
        self.sge = sge
        self.sge_opts = sge_opts
        
        # Setup logging
        if not op.isdir(log_dir):
            os.mkdir(log_dir)
        dt = datetime.now()
        self.logdir = op.join(log_dir, dt.strftime("%Y-%m-%d_%H:%M"))
        if not op.isdir(self.logdir):
            os.mkdir(self.logdir)
        log = op.join(self.logdir, "niwrapper.log")
        ## save command that run for this
        fname = op.join(self.logdir, "orig_command.txt")
        f = file(fname, 'w')
        f.write(" ".join(self.argv))
        f.close()
        ## save subjects
        fname = op.join(self.logdir, "subjects.txt")
        f = file(fname, 'w')
        f.write("\n".join(subjects))
        f.close()
        
        # Call parent
        super(NiWrapper, self).__init__(verbosity, self.template_context, dry_run, log)
        self.log.debug("NiWrapper logger is up and running")
        
        self._is_parsed = True
        return
    
    def parse_parser(self, arglist):
        self._create_parser()
        args = self.parser.parse_args(arglist)
        self.setup(**args)
        return
    
    def _create_parser(self, **kwargs):
        parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS, **kwargs)
        
        group = parser.add_argument_group("Command Options")
        for k,v in self._parser_help.iteritems():
            group.add_argument("--%s" % k, action="append_const", const=k, dest="run_keys", help=v)
        
        group = parser.add_argument_group('Subject Options')
        group.add_argument('-s', '--subjects', nargs="+", required=True, metavar="subject")
        group.add_argument('-q', '--sge', action="store_true", default=False, help="use SGE")
        group.add_argument('--sge-opts', default="-S /bin/bash -V -cwd -j y", 
                            help="default: %(default)s")
        
        group = parser.add_argument_group('I/O Options')
        group.add_argument("--verbose", action="store_const", const=1, dest="verbosity", default=0)
        group.add_argument("--debug", action="store_const", const=2, dest="verbosity", default=0)
        group.add_argument("--dry-run", action="store_true", default=False)
        group.add_argument('--log-dir', default=op.join(os.getcwd(), "nilogs"),
                            help="default: %(default)s")
        
        self.parser = parser
        return
    
    def compile(self):
        if not self._is_parsed:
            raise Exception("Have not parsed anything yet")
        
        self.log.info("Compiling")
        commands = {}
        for k,v in self._commands_opts.iteritems():
            cmd = []
            prog,opts = v
            for ko,vo in opts:
                if len(ko) == 1:
                    pre = "-"
                else:
                    pre = "--"
                if isinstance(vo, bool) and vo == True:
                    cmd.append("%s%s" % (pre, ko))
                else:
                    cmd.append("%s%s %s" % (pre, ko, self._substitute(vo)))
            commands[k] = "%s %s" % (prog, " ".join(cmd))
        self._commands = commands
        return commands
    
    def run(self):
        if not self._is_parsed:
            raise Exception("Have not parsed anything yet")
        
        self.compile()
        self.log.info("Running")
        
        self.log.debug("checking command labels")
        for k in self.run_keys:
            if k not in self._commands:
                self.log.error("Could not find command label %s" % k)
        
        self.log.debug("loop through each participant and run command")
        if self.sge:
            self.log.debug("...using SGE")
            self._setup_sge()
        for s in subject_list:
            self.log.title("Subject: %s" % s)
            for k in self.run_keys:
                self.log.subtitle("command: %s" % k)
                cmd = "%s -s %s" % (self._commands[k], s)
                self._execute(s, k, cmd)
        
        return
    
    def _setup_sge(self):
        self.sge_stdout = op.join(self.logdir, "sge_outputs")
        self.sge_scripts = op.join(self.logdir, "sge_scripts")
        if not op.isdir(self.sge_stdout):
            os.mkdir(self.sge_stdout)
        if not op.isdir(self.sge_scripts):
            os.mkdir(self.sge_scripts)
        return
    
    def _execute(self, subject, label, cmd):
        if not self._is_parsed:
            raise Exception("Have not parsed anything yet")
        if self.sge:
            script = op.join(self.sge_scripts, "%s_%s.bash")
            output = op.join(self.sge_stdout, "%s_%s.txt")
            # Create a new file
            f = file(script, 'w')
            f.write("#!/usr/bin/env bash\n")
            f.write("\n")
            f.write(cmd + "\n")
            f.write("\n")
            f.close()
            # New command
            cmd = "qsub %s -o %s %s" % (self.sge_opts, output, script)
        # Execute
        p = Popen(cmd, shell=True, cwd=os.getcwd())
        p.communicate()
        if p.returncode != 0:
            self.log.error("Error running '%s'" % cmd)
        return p.returncode
    

def main(arglist):
    # Find config file in current directory
    global config_file
    if len(arglist) > 0 and arglist[0][0] != '-' and arglist[0][0] != '--':
        config_file = arglist[0]
        arglist = arglist[1:]
        print 'Manually defined config file: %s' % config_file
    if not op.isfile(config_file):
        print("ERROR: config file '%s' was not found in current directory" % config_file)
        raise SystemExit(2)
    niwrap = NiWrapper(config_file)
    niwrap.parse_parser(arglist)
    niwrap.run()

if __name__ == "__main__":
    main(sys.argv[1:])
