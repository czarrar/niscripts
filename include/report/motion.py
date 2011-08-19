import os, sys
import os.path as op
from glob import glob
from jinja2 import Environment, PackageLoader

sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))
import zlogger

class MotionReporter(object):
    """docstring for MotionReporter"""
    _logname = "motion_report"
    _motion_max_fname = "motion_max.txt"
    _motion_disp_fname = "motion_disp.png"
    
    def __init__(self, verbosity, log=None):
        super(MotionReport, self).__init__()
        
        if verbosity == 0:
            self.loglevel = zlogger.logging.IMPORTANT
        elif verbosity == 1:
            self.loglevel = zlogger.logging.COMMAND_INFO
        else:
            self.loglevel = zlogger.logging.DEBUG
        self.log = zlogger.getLogger(self._logname, level=self.loglevel, fname=log, 
                                     allow_exceptions=True)
        self.log.info("Starting Motion Reporter")
        
        self.env = Environment(loader=PackageLoader('report', 'templates'))
        self.subjects = {}
        self.report = None
        
        return
    
    def compile(self):
        self.log.info("Compiling Report Page")
        if not self.subjects:
            self.log.fatal("Subjects not set!")
        if report:
            self.log.warn("Report was already set")
        context = {'title': 'Motion', 'subjects': self.subjects}
        self.report = self._render('motion', **context)
        return self.report
    
    def run(self, redo=False):
        if not self.report and not redo:
            self.compile()
        self.log.info("Writing %s" % self.outfile)
        f = file(self.outfile, 'w')
        f.write(self.report)
        f.close()
        return
    
    def _remove_extra_stuff(self, a):
        b = re.sub("\n\ +", "\n", a)
        c = re.sub("\n\n+", "\n\n", b)
        return c
    
    def _render(self, fname, **context):
        template = self.env.get_template(fname + '.jinja')
        result = template.render(**context)
        return self._remove_extra_stuff(result)
    
    def addSubject(self, subject, indirs):
        self.log.debug("...adding subject %s" % subject)
        sinfo = {}
        for run,indir in enumerate(indirs):
            if not op.path(indir):
                self.log.error("Input directory '%s' does not exist" % indir)
            # new dict
            rinfo = {}
            # max motion
            fmax = op.join(indir, self._motion_max_fname)
            if not op.isfile(fmax):
                self.log.error("Max motion file '%s' does not exist" % fmax)
            f = file(infile, 'r')
            lines = f.readlines()
            f.close()
            rinfo['abs'] = lines[1].strip()
            rinfo['rel'] = lines[3].strip()
            # disp pic
            fdisp = op.join(indir, self._motion_disp_fname)
            if not op.isfile(fdisp):
                self.log.error("Displacment picture '%s' does not exist" % fdisp)
            rinfo['disp'] = fdisp
            # save
            sinfo[run] = rinfo
        self.subjects[subject] = sinfo        
        return
    
    def setData(self, subject_indirs, outfile):
        self.log.info("Setting Data")
        for subject,indir in subject_indirs.iteritems():
            indirs = glob(op.expanduser(indir))
            if len(indirs) == 0:
                self.log.error("Couldn't find '%s' for %s" % (indir, subject))
            self.addSubject(subject, indirs)
        self.outfile = op.splitext(outfile)[0] + ".html"
        return
    

