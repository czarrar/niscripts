#!/usr/bin/env python

import argparse, os, yaml, sys
import os.path as op
from datetime import datetime
from analysis.base import SubjectBase
from subprocess import Popen
#from zlogger import (LoggerError, LoggerCritical)
from usage import append_var
from glob import glob
from copy import deepcopy

config_file = "niwrapper.yaml"

def die(msg):
    print "ERROR: %s" % msg
    raise SystemExit(2)

def tolist(x):
    if not isinstance(x, list):
        x = [x]
    return x

class NiWrapper(SubjectBase):
    # Scripts that make use of nipype
    _use_nipype = [
        "preprocess_anat.py", "preprocess_func.py", "register.py", "nuisance_evs.py", 
        "preprocess_func_p2.py", 
    ]
    # Scripts in niscripts (that use the -s option)
    _niscripts = [ op.basename(x) 
                    for x in glob(op.join(os.getenv("NISCRIPTS"),"commands","*.py")) ]
    
    def __init__(self, config_file):
        # Load config
        if not op.isfile(config_file):
            raise Exception("Cannot find configuration file: %s" % config_file)
        f = file(config_file, 'r')
        self.config = yaml.load(f)
        f.close()
        
        # Get context
        self.template_context = self.config.pop("vars", {})
        
        # Check that certain variables are given via the command-line
        self.require_user_vars = self.config.pop("require_user_vars", [])
        
        # Setup
        self._setup_command_opts()
        
        self._is_parsed = False
        return
    
    def _setup_command_opts(self):
        self._parser_help = {}
        self._commands_opts = {}
        self._short_opt = {}
        self._long_opt = {}        
        self._inputs = {}
        self._outputs = {}
        self._overwrite = {}
        self._create_dirs = {}
        for k,opts in self.config.iteritems():
            try:
                params = opts.pop("_params")
            except KeyError:
                die("You must specify the program paramaters (_params) in %s" % k)
            # help
            if k in self._parser_help:
                die("Duplicate sub-command %s" % k)
            self._parser_help[k] = params.pop("help", "sub-command")
            # command
            if k in self._commands_opts:
                die("Duplicate sub-command %s" % k)
            try:
                p = params.pop("program")
            except KeyError:
                die("Must specify program in config file (within _params of %s)" % k)
            # short_opt/long_opt (e.g., - or -- or +)
            self._short_opt[k] = params.pop('short_opt', '-')
            self._long_opt[k] = params.pop('long_opt', '--')
            # inputs to check
            self._inputs[k] = {}
            if 'in_opts' in params:
                inopts = tolist(params.pop('in_opts'))
                for inopt in inopts:
                    if inopt not in opts:
                        die("%s specified in in_opts but not an option in %s" % 
                                (inopt, k))
                    self._inputs[k].extend(str(opts[inopt]))
            if 'inputs' in params:
                self._inputs[k].append(tolist(params.pop('inputs')))
            # outputs to check
            self._outputs[k] = {}
            if 'out_opts' in params:
                outopts = tolist(params.pop('out_opts'))
                for outopt in outopts:
                    if outopt not in opts:
                        die("%s specified in out_opts but not an option in %s" % 
                                (outopt, k))
                    self._outputs[k].append(str(opts[outopt]))
            if 'outputs' in params:
                self._outputs[k].extend(tolist(params.pop('outputs')))
            # overwrite?
            self._overwrite[k] = params.pop('overwrite', False)
            # ensure certain directories exist before running program
            self._create_dirs[k] = tolist(params.pop('create_dirs', []))
            # check params
            if params:
                die('some stuff is still left in params for %s: %s' % (k, params))
            # save options
            self._commands_opts[k] = (p, opts)
        return self._commands_opts
    
    def setup(self, run_keys, log_dir, subjects, sge, sge_opts, verbosity, processors, 
                dry_run=False, vars={}, **cmds):
        # Set subjects
        new_subjects = []
        for s in subjects:
            if s[0] == "@":
                f = file(s[1:], 'r')
                tmps = f.readlines()
                tmps = [ l.strip() for l in tmps if l.strip() and l.strip()[0] != "#" ]
                new_subjects.extend(tmps)
            else:
                new_subjects.append(s)
        subjects = new_subjects
        
        # Save
        self.template_context.update(vars)
        for x in self.require_user_vars:
            if x not in self.template_context:
                die("The variable %s is required but not specified (use --var %s=XXX)" % 
                    (x, x))
        self.processors = processors
        self.run_keys = run_keys
        self.subjects = subjects
        self.sge = sge
        self.sge_opts = sge_opts
        
        # Setup logging
        if not op.isdir(log_dir):
            os.mkdir(log_dir)
        dt = datetime.now()
        self.logdir = op.join(log_dir, dt.strftime("%Y-%m-%d_%H-%M") + "_" + "_".join(run_keys))
        if not op.isdir(self.logdir):
            os.mkdir(self.logdir)
        log = op.join(self.logdir, "niwrapper.log")
        ## save command that run for this
        fname = op.join(self.logdir, "orig_command.txt")
        f = file(fname, 'w')
        f.write(" ".join(sys.argv))
        f.close()
        ## save subjects
        fname = op.join(self.logdir, "subjects.txt")
        f = file(fname, 'w')
        f.write("\n".join(subjects))
        f.close()
        
        # Call parent
        super(NiWrapper, self).__init__(verbosity, self.template_context, dry_run, log)
        self.log.debug("NiWrapper logger is up and running")
        
        # Check
        if self.sge and self.processors > 1:
            self.log.error("Cannot specify -q/--sge and -p/--processors")
        
        self._is_parsed = True
        return
    
    def parse_parser(self, arglist):
        self._create_parser()
        args = self.parser.parse_args(arglist)
        self.setup(**vars(args))
        return
    
    def _create_parser(self, **kwargs):
        parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS, **kwargs)
        
        group = parser.add_argument_group("Command Options")
        for k,v in self._parser_help.iteritems():
            group.add_argument("--%s" % k, action="append_const", const=k, dest="run_keys", help=v)
        
        group = parser.add_argument_group('Subject Options')
        group.add_argument('-s', '--subjects', nargs="+", required=True, metavar="subject")
        group.add_argument('--var', action="append", type=append_var, dest="vars")
        group.add_argument('-p', '--processors', default=1, type=int, help="use multiple processers on same computer (must specify number of processors)")
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
        self._workingdirs = {}
        commands_opts = deepcopy(self._commands_opts)
        for k,v in commands_opts.iteritems():
            cmd = []
            prog,opts = v
            if self.sge:
                w = opts.pop('workingdir', '')
                self._workingdirs[k] = self._substitute(w)
            for ko,vo in opts.iteritems():
                if len(ko) == 1:
                    pre = self._short_opt[k]
                else:
                    pre = self._long_opt[k]
                if isinstance(vo, bool) and vo == True:
                    cmd.append("%s%s" % (pre, ko))
                elif isinstance(vo, list):
                    for ivo in vo:
                        cmd.append("%s%s %s" % (pre, ko, self._substitute(str(ivo))))
                else:
                    cmd.append("%s%s %s" % (pre, ko, self._substitute(str(vo))))
            if prog in self._use_nipype:
                cmd.append("--crash-dir %s" % op.join(self.logdir, "crashes"))
            commands[k] = "%s %s" % (prog, " ".join(cmd))
        self._commands = commands
        return commands
    
    def run(self):
        if not self._is_parsed:
            raise Exception("Have not parsed anything yet")
        
        self.log.info("Running")
                
        self.log.debug("loop through each participant and run command")
        if self.sge:
            self.log.debug("...using SGE")
            self._setup_sge()
        if self.processors > 1:
            self.compile()
            self.log.debug("checking command labels")
            for k in self.run_keys:
                if k not in self._commands:
                    self.log.error("Could not find command label %s" % k)
            for k in self.run_keys:
                self.log.subtitle("command: %s" % k)
                cmd = "%s --plugin MultiProc %i -s %s" % (self._commands[k], self.processors, 
                                                            " ".join(self.subjects))
                self._execute(cmd, label=k)
        else:
            for s in self.subjects:
                self.log.title("Subject: %s" % s)
                self.template_context['subject'] = s
                self.compile()
                self.log.debug("checking command labels")
                for k in self.run_keys:
                    if k not in self._commands:
                        self.log.error("Could not find command label %s" % k)
                for k in self.run_keys:
                    self.log.subtitle("command: %s" % k)
                    self._execute(self._commands[k], s, k)
                del self.template_context['subject']
        return
    
    def _setup_sge(self):
        self.sge_stdout = op.join(self.logdir, "sge_outputs")
        self.sge_scripts = op.join(self.logdir, "sge_scripts")
        if not op.isdir(self.sge_stdout):
            os.mkdir(self.sge_stdout)
        if not op.isdir(self.sge_scripts):
            os.mkdir(self.sge_scripts)
        return
    
    def _execute(self, cmd_opt, subject=None, label=None):
        if not self._is_parsed:
            raise Exception("Have not parsed anything yet")
        prog = cmd_opt.split(" ")[0]
        if self._create_dirs[label]:
            for d in self._create_dirs[label]:
                if not op.isdir(d):
                    self.log.debug("creating directory %s" % d)
                    os.mkdir(d)
        if self.sge:
            if not self._workingdirs[label]:
                print self._workingdirs
                if prog in self._use_nipype:
                    self.log.error("You must specify the workingdir for %s" % prog)
            else:
                wd = op.join(self._workingdirs[label], subject)
                if not op.isdir(wd):
                    os.mkdir(wd)
                cmd_opt = "%s --workingdir %s" % (cmd_opt, wd)
            if prog in self._niscripts:
                cmd = "%s -s %s" % (cmd_opt, subject)
            if subject is None or label is None:
                self.log.fatal("Must specificy subject and label for _execute")
            script = op.join(self.sge_scripts, "x_%s_%s.bash" % (subject, label))
            output = op.join(self.sge_stdout, "o_%s_%s.txt" % (subject, label))
            # Create a new file
            f = file(script, 'w')
            f.write("#!/usr/bin/env bash\n")
            f.write("\n")
            f.write(cmd + "\n")
            f.write("\n")
            f.close()
            # New command
            cmd = "qsub %s -o '%s' '%s'" % (self.sge_opts, output, script)
        elif subject is not None:
            if prog in self._niscripts:
                cmd = "%s -s %s" % (cmd_opt, subject)
            else:
                cmd = "%s" % cmd_opt
        else:
            cmd = "%s" % cmd_opt
        # Execute
        self.log.drycommand(cmd)
        if not self.dry_run:
            p = Popen(cmd, shell=True, cwd=os.getcwd(), stdout=sys.stdout, stderr=sys.stderr)
            p.communicate()
            if p.returncode != 0:
                self.log.error("Error running command")
            return p.returncode
        else:
            return 0
    

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
