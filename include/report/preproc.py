import os, re, sys
import os.path as op
from glob import glob
from jinja2 import Environment, PackageLoader

sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))
import zlogger
from report.base import Reporter, ReporterException

class PreprocReporter(Reporter):
    """docstring for PreprocReporter"""
    _logname = "preproc_report"
    _motion_max_fname = "motion_max.txt"
    _motion_disp_fname = "motion_disp.png"
    
    def __init__(self, *args, **kwargs):
        super(PreprocReporter, self).__init__(*args, **kwargs)
        self.log.title("Starting Preprocessing Reporter")        
        self.subjects = []
        return
    
    def run(self):
        if not self.subjects:
            self.log.fatal("Subjects not set!")
                
        self.log.info("Creating report")
        for subject,context in self.subjects:
            self.log.info("...subject %s" % subject)
            subfile = op.join(self.outdir, subject + ".html")
            context['title'] = 'Preprocessing'
            report = self._render('preproc_subjects', **context)        
            self.log.debug("...writing %s" % subfile)
            f = file(subfile, 'w')
            f.write(report)
            f.close()
        return
    
    def _psub(self, path):
        return op.relpath(self._substitute(x), self.outdir)
    
    def addSubject(self, subject, anatdirs, regdirs, funcdirs):
        self.log.info("...subject: %s" % subject)
        if len(regdirs) != len(funcdirs):
            self.log.fatal("Number of functional and run directories must be the same")
        if regdirs.keys() != funcdirs.keys():
            self.log.fatal("Functional directory and run directory labels must be the same")
        
        sinfo = {}
        ainfo = {}  # mask
        rinfo = {}  # reg
        minfo = {}  # motion
        
        # setup
        self.log.debug("...setup")
        self.addTemplateContext('subject', subject)
        ## anat
        anatdirs = [ (l,self._psub(x)) for l,x in anatdirs.iteritems() ]
        [ 
            self.log.error("anatomical directory '%s' does not exist" % x)
            for l,x in anatdirs if not op.isdir(x)
        ]
        ## reg
        regdirs = [ (l,self._psub(x)) for l,x in regdirs.iteritems() ]
        [ 
            self.log.error("registration directory '%s' does not exist" % x)
            for l,x in regdirs if not op.isdir(x)
        ]
        ## func
        tmp_funcdirs = [ (l,glob(self._psub(x))) for l,x in funcdirs.iteritems() ]
        [ 
            self.log.error("functional directory '%s' does not exist" % funcdirs[l])
            for l,x in tmp_funcdirs if len(x)==0
        ]
        funcdirs = tmp_funcdirs
        ## save
        sinfo['anatdirs'] = anatdirs
        sinfo['regdirs'] = regdirs
        sinfo['funcdirs'] = funcdirs
        
        # anat
        self.log.debug("...anat")
        for label,anatdir in anatdirs:
            ainfo[label] = [
                ('head_axial', op.join(anatdir, "head_axial.png")), 
                ('mask_axial', op.join(anatdir, "brain_mask_axial.png")), 
                ('mask_sagittal', op.join(anatdir, "brain_mask_sagittal.png"))
            ]
            for k,v in ainfo[label]:
                if not op.isfile(v):
                    self.log.error("anatomical pic '%s' does not exist" % v)
        
        # reg
        self.log.debug("...reg")
        for label,regdir in regdirs:
            rinfo[label] = []
            if op.isfile(op.join(regdir, "func2coplanar.png")):
                rinfo[label].append(("func2coplanar", op.join(regdir, "func2coplanar.png")))
            if op.isfile(op.join(regdir, "coplanar2highres.png")):
                rinfo[label].append(("coplanar2highres", op.join(regdir, "coplanar2highres.png")))
            rinfo[label].extend([
                ("func2highres", op.join(regdir, "func2highres.png")), 
                ("highres2standard_lin", op.join(regdir, "highres2standard.png")), 
                ("func2standard_lin", op.join(regdir, "func2standard.png"))
            ])
            if op.isfile(op.join(regdir, "highres2standard_fnirt.png")):
                rinfo[label].append(("highres2standard_nonlin", 
                                        op.join(regdir, "highres2standard_fnirt.png")))
            if op.isfile(op.join(regdir, "func2standard_fnirt.png")):
                rinfo[label].append(("func2standard_nonlin", 
                                        op.join(regdir, "func2standard_fnirt.png")))
            for k,v in rinfo[label]:
                if not op.isfile(v):
                    self.log.error("registration pic '%s' does not exist" % v)
        
        # motion
        self.log.debug("...motion")
        for label,rundirs in funcdirs:
            minfo[label] = []
            for run,rundir in enumerate(rundirs):
                # new dict
                tmpinfo = {}
                # max motion
                fmax = op.join(rundir, self._motion_max_fname)
                if not op.isfile(fmax):
                    self.log.error("Max motion file '%s' does not exist" % fmax)
                f = file(fmax, 'r')
                lines = f.readlines()
                f.close()
                tmpinfo['abs'] = "%.2f" % float(lines[1].strip())
                tmpinfo['rel'] = "%.2f" % float(lines[3].strip())
                # disp pic
                fdisp = op.join(rundir, self._motion_disp_fname)
                if not op.isfile(fdisp):
                    self.log.error("Displacment picture '%s' does not exist" % fdisp)
                tmpinfo['disp'] = fdisp
                # save
                minfo[label].append((run+1,tmpinfo))
        
        self.log.debug("...saving")
        self.addTemplateContext('subject', None)
        self.subjects.append((subject, {
            'anat': ainfo,
            'reg': rinfo,
            'motion': minfo
        }))
        return
    
    def setData(self, subjects, anatdirs, regdirs, funcdirs):
        self.log.info("Setting Data")
        self.subjects = []
        for subject in subjects:
            self.addSubject(subject, anatdirs, regdirs, funcdirs)
        return
    

