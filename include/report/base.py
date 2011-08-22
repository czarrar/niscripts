import os, shutil, sys
import os.path as op

sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))
import zlogger

class ReporterException(Exception):
    """todo"""

class Reporter(object):
    """docstring for Reporter"""
    _logname = None
    
    def __init__(self, verbosity, outdir, template_context={}, overwrite=False):
        super(PreprocReporter, self).__init__()
        
        outdir = op.abspath(op.expanduser(outdir))
        if op.isdir(outdir):
            if overwrite:
                shutil.rmtree(outdir)
            else:
                raise ReporterException("Output directory already exists")
        os.mkdir(outdir)
        self.outdir = outdir
        os.mkdir(op.join(outdir, "logs"))
        if self._logname is None:
            raise Exception("You must set the _logname in the class")
        logfile = op.join(op.join(outdir, "logs"), self._logname + ".log")
        
        if verbosity == 0:
            self.loglevel = zlogger.logging.IMPORTANT
        elif verbosity == 1:
            self.loglevel = zlogger.logging.COMMAND_INFO
        else:
            self.loglevel = zlogger.logging.DEBUG
        self.log = zlogger.getLogger(self._logname, level=self.loglevel, fname=logfile, 
                                     allow_exceptions=True)
        
        self._setTemplateContext(template_context)
        self.env = Environment(loader=PackageLoader('report', 'templates'))
        self.report = None
        
        return
    
    def addTemplateContext(self, k, v):
        self.template_context[k] = v
    
    def _setTemplateContext(self, context):
        if not context:
            self.template_context = {}
            return
        for i in xrange(10):
            to_stop = True
            for k,v in context.iteritems():
                if not v:
                    self.log.fatal("Setting template: value for variable '%s' is empty" % k)
                t = Template(v)
                try:
                    context[k] = t.substitute(**context)
                    if (len(list_template_vars(v)) > 0):
                        to_stop = False
                except KeyError as kk:
                    self.log.fatal("Setting template: could not find value for '%s' in '%s'" 
                                        % (kk,k))
                    context[k] = t.safe_substitute(**context)
            if to_stop:
                break
        if not to_stop:
            self.log.warning("Setting template: did not finish substitutions of template context")
        self.template_context = context
        return
    
    def _substitute(self, template, maxloops=5):
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
    
    def _remove_extra_stuff(self, a):
        b = re.sub("\n\ +", "\n", a)
        c = re.sub("\n\n+", "\n\n", b)
        return c
    
    def _render(self, fname, **context):
        template = self.env.get_template(fname + '.jinja')
        result = template.render(**context)
        return self._remove_extra_stuff(result)
    

