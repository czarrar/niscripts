import zlogger, re, os, shutil
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

def get_loglevel(verbosity):
    if verbosity == 0:
        loglevel = zlogger.logging.IMPORTANT
    elif verbosity == 1:
        loglevel = zlogger.logging.COMMAND_INFO
    else:
        loglevel = zlogger.logging.DEBUG
    return loglevel

def create_logger(name, level, fname):
    log = zlogger.getLogger(name, level, fname=fname, allow_exceptions=True)
    return log

class Base(object):
    _logname = "default"
    
    def __init__(self, verbosity, template_context, dry_run=False, log=None, logger=None):
        """
        Parameters
        ----------
        log_name
        """
        super(Base, self).__init__()
        
        self.loglevel = get_loglevel(verbosity)
        
        # ghetto fix!!!
        if not logger:
            self.log = create_logger(self._logname, self.loglevel, log)
        else:
            self.log = logger
        
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
    
    def _input_list2str(self, inputs, desc="input files"):
        if len(inputs) > 1:
            self.log.error("Too many %s found: %s" % (decs, infunc))
        return inputs[0]
    
    def _check_infile(self, fpath, desc="input file", substitute=False):
        if substitute:
            fpath = self._substitute(fpath)
        if not op.isfile(fpath):
            self.log.error("%s '%s' does not exist." % (desc, fpath))
        return fpath
    
    def _check_indir(self, dpath, desc="input directory", substitute=False):
        if substitute:
            dpath = self._substitute(dpath)
        if not op.isdir(dpath):
            self.log.error("%s '%s' does not exist." % (desc, dpath))
        return dpath
    
    def _check_outfile(self, fpath, desc="output file", overwrite=False, substitute=False):
        if substitute:
            fpath = self._substitute(fpath)
        if op.isfile(fpath):
            if overwrite:
                self.log.warning("removing %s '%s'" % (desc, fpath))
                os.remove(fpath)
            else:
                self.log.error("%s '%s' already exists" % (desc, fpath))
        return fpath
    
    def _check_outdir(self, dpath, desc="output directory", overwrite=False, substitute=False):
        if substitute:
            dpath = self._substitute(dpath)
        if op.isfile(dpath):
            if overwrite:
                self.log.warning("removing %s '%s'" % (desc, dpath))
                shutil.rmtree(dpath)
            else:
                self.log.error("%s '%s' already exists" % (desc, dpath))
        return dpath
    
    def _getInputsWorker(self, infiles, itype='path'):
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
           new_infiles.sort()
           return new_infiles
    
    def addTemplateContext(self, k, v):
        self.template_context[k] = v
    
    def delTemplateContext(self, k):
        del self.template_context[k]
    
    def _setTemplateContext(self, context):
        self.log.debug("setting template context")
        for k,v in context.iteritems():
            self.log.debug("%s = %s" % (k,v))
            if not v:
                self.log.fatal("Value for variable '%s' is empty" % k)
            t = Template(v)
            try:
                context[k] = t.substitute(**context)
            except KeyError as kk:
                self.log.debug("Could not find value for '%s' in '%s'" % (kk,k))
                context[k] = t.safe_substitute(**context)
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
    


class SubjectBase(Base):
    """
    Base Class for CombineFuncs and FsfSubject
    
    Examples
    --------
    >>> from feat_subject import SubjectBase
    >>> sb = SubjectBase('test', 2, {'base': "/path/to/base", 'jack': "${base}/jack", 
                                    'black': "${jack}/black"})
    >>> sb._substitute("${black}")   # should have /path/to/base/jack/black
    """
    
    def _file_ncols(self, fname):
        f = open(fname, 'r')
        l = f.readline()
        f.close()
        ncols = len(re.split("[ \t]+", l.strip()))
        return ncols
    
    def _file_nrows(self, fname):
        f = open(fname, 'r')
        ls = f.readlines()
        f.close()
        nrows = len([ l for l in ls if l.strip() ])
        return nrows
    
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
    

