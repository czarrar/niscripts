import zlogger, re, os
import os.path as op
from glob import glob
from string import Template

import sys
sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))
from execute import Process

def list_jinja2_vars(fprefix):
    """Utility function for testing purposes"""
    full_fname = op.join(os.getenv("NISCRIPTS"), "include", "feat", "templates", 
                            fprefix + ".jinja")
    f = file(full_fname, 'r')
    s = f.read()
    return list(set(re.findall("{{\ *(.*?)\ *}}", s)))

def list_template_vars(string):
    template = Template(string)
    tmp = template.pattern.findall(template.template)
    res = []
    for x in tmp:
        if x[1]:
            res.append(x[1])
        elif x[2]:
            res.append(x[2])
    return res


class SubjectBase(object):
    """
    Base Class for CombineFuncs and FsfSubject
    
    Examples
    --------
    >>> from feat_subject import SubjectBase
    >>> sb = SubjectBase('test', 2, {'base': "/path/to/base", 'jack': "${base}/jack", 
                                    'black': "${jack}/black"})
    >>> sb._substitute("${black}")   # should have /path/to/base/jack/black
    """
    
    _logname = "default"
    
    def __init__(self, verbosity, template_context, dry_run=False, log=None):
        """
        Parameters
        ----------
        log_name
        """
        super(SubjectBase, self).__init__()
        
        if verbosity == 0:
            self.loglevel = zlogger.logging.IMPORTANT
        elif verbosity == 1:
            self.loglevel = zlogger.logging.COMMAND_INFO
        else:
            self.loglevel = zlogger.logging.DEBUG
        self.log = zlogger.getLogger(self._logname, level=self.loglevel, fname=log, allow_exceptions=True)
        
        self._setTemplateContext(template_context)
        
        self.dry_run = dry_run
        
        return
    
    def check_req(self, d, l):
        """Check that each element of l is in d"""
        for e in l:
            if e not in d:
                self.log.error("%s must be in options" % e)

    def check_ex(self, d, l):
        """Check that each element of l is NOT in d"""
        for e in l:
            if e in d:
                self.log.error("%s must not be in options" % e)
    
    def _file_ncols(self, fname):
        f = open(fname, 'r')
        l = f.readline()
        f.close()
        ncols = len(re.split("[ \t]+", l.strip()))
        return ncols
    
    def addTemplateContext(self, k, v):
        self.template_context[k] = v
    
    def _setTemplateContext(self, context):
        self.log.debug("setting template context")
        for i in xrange(10):
            to_stop = True
            for k,v in context.iteritems():
                if not v:
                    self.log.fatal("Value for variable '%s' is empty" % k)
                t = Template(v)
                try:
                    context[k] = t.substitute(**context)
                    if (len(list_template_vars(v)) > 0):
                        to_stop = False
                except KeyError as kk:
                    self.log.warning("Could not find value for '%s' in '%s'" % (kk,k))
                    context[k] = t.safe_substitute(**context)
            if to_stop:
                break
        if not to_stop:
            self.log.warning("Did not finish substitutions of template context")
        self.template_context = context
    
    def _substitute(self, template, maxloops=5):
        self.log.debug("substituting %s" % template)
        leftover = 0
        result = template
        for i in xrange(maxloops):
            try:
                result = Template(result).substitute(**self.template_context)
            except KeyError as k:
                self.log.fatal("Could not find value for %s in '%s'" % (k, template))
            leftover = len(list_template_vars(result))
            if (leftover == 0):
                break
        if leftover > 0:
            self.log.warning("Did not finish substitition of '%s'" % template)
        return result
    
    def _getInputsWorker(self, infiles, itype):
        new_infiles = []
        if isinstance(infiles, str):
            infiles = self._substitute(infiles)
            new_infiles = glob(infiles)
            if len(new_infiles) == 0:
                self.log.error("Input %s '%s' does not exist" % (itype, infiles))
        elif isinstance(infiles, list):
            for infile in infiles:
                new_infile = glob(self._substitute(infile))
                if len(new_infile) == 0:
                    self.log.error("List input %s '%s' does not exist" % (itype, infile))
                elif len(new_infile) > 1:
                    self.log.warning("List input '%s' points to more than one %s" % (infile, 
                                                                                        itype))
                new_infiles.extend(new_infile)
        else:
            raise Exception("Incorrect type for infiles: %s" % infiles)
        return new_infiles
    
    def getInputs(self, infiles, runfile, mkdir):
        # Run file
        if runfile:
            runfile = self._substitute(runfile)
            if not op.isfile(runfile):
                self.log.error("Could not find run file '%s'" % runfile)
            f = file(runfile, 'r')
            runs = f.readlines()
            runs = [ l.strip() for l in runs if l.strip() ]
            if len(runs) == 0:
                self.log.error("Empty run file '%s'" % runfile)
        else:
            runs = None
        
        # Set infiles
        new_infiles = []
        if runs:
            for r in runs:
                self.addTemplateContext('run', r)
                new_infiles.extend(self._getInputsWorker(infiles, 'file'))
            self.addTemplateContext('run', None)
        else:
            new_infiles.extend(self._getInputsWorker(infiles, 'file'))
        
        # Make a directory
        if mkdir:
            if isinstance(mkdir, str):
                mkdir = [mkdir]
            for m in mkdir:
                m = self._substitute(m)
                if m and not op.isdir(m):
                    self.log.info("Making directory: %s" % m)
                    os.mkdir(m)
        
        self.runs = runs
        return new_infiles
    
    def setMotion(self, infiles, outfile, mkdir=None):
        # Set infiles
        new_infiles = []
        if self.runs:
            for r in self.runs:
                self.addTemplateContext('run', r)
                new_infiles.extend(self._getInputsWorker(infiles, 'file'))
            self.addTemplateContext('run', None)
        else:
            new_infiles.extend(self._getInputsWorker(infiles, 'file'))
        
        # Make a directory
        if mkdir:
            if isinstance(mkdir, str):
                mkdir = [mkdir]
            for m in mkdir:
                m = self._substitute(m)
                if m and not op.isdir(m):
                    self.log.info("Making directory: %s" % m)
                    os.mkdir(m)
        
        outfile = self._substitute(outfile)
        self.log.info("Creating new combined motion file '%s'" % outfile)
        f = file(outfile, 'w')
        p = Process("cat %s" % " ".join(new_infiles), stdout=f, to_print=True)
        print p.stderr
        self.outmotion = outfile
        
        return
    

