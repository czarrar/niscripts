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
    _motion_outliers_fname = "pics_outlier_vols.png"
    
    def __init__(self, *args, **kwargs):
        super(PreprocReporter, self).__init__(*args, **kwargs)
        self.log.title("Starting Preprocessing Reporter")        
        self.subjects = []
        return
    
    def run(self):
        if not self.subjects:
            self.log.fatal("Subjects not set!")
        
        self.log.info("Creating report")
        links = [ (s, s + ".html") for s,c in self.subjects ]
        for subject,context in self.subjects:
            self.log.info("...subject %s" % subject)
            subfile = op.join(self.outdir, subject + ".html")
            context['subject'] = subject
            context['title'] = 'Preprocessing'
            context['links'] = links
            report = self._render('preproc_subjects', **context)        
            self.log.debug("...writing %s" % subfile)
            f = file(subfile, 'w')
            f.write(report)
            f.close()
        
        # index page
        self.log.info("...index page")
        ofile = op.join(self.outdir, "index.html")
        context['title'] = 'Preprocessing'
        context['subjects'] = links
        report = self._render('preproc_index', **context)
        self.log.debug("...writing %s" % ofile)
        f = file(ofile, 'w')
        f.write(report)
        f.close()
        
        return
    
    def _psub(self, path):
        return op.relpath(self._substitute(path), self.outdir)
    
    def addSubject(self, subject, anatdirs, funcdirs, regdirs, segdirs, nuidirs):
        cwd = os.getcwd()
        os.chdir(self.outdir)
        self.log.info("...subject: %s" % subject)
        
        sinfo = {}
        ainfo = {}  # mask
        sinfo = {}  # segmentation
        minfo = {}  # motion
        rinfo = {}  # reg
        ninfo = {}  # nuisance
        
        # setup
        self.log.debug("...setup")
        self.addTemplateContext('subject', subject)
        ## anat
        anatdirs = [ (l,self._psub(x)) for l,x in anatdirs.iteritems() ]
        [ 
            self.log.error("anatomical directory '%s' does not exist" % x)
            for l,x in anatdirs if not op.isdir(x)
        ]
        ## segmentation
        if segdirs:
            segdirs = [ (l,self._psub(x)) for l,x in segdirs.iteritems() ]
            [ 
                self.log.error("segmentation directory '%s' does not exist" % x)
                for l,x in segdirs if not op.isdir(x)
            ]
        ## func
        if funcdirs:
            tmp_funcdirs = [ (l,glob(self._psub(x))) for l,x in funcdirs.iteritems() ]
            [ 
                self.log.error("functional directory '%s' does not exist" % funcdirs[l])
                for l,x in tmp_funcdirs if len(x)==0
            ]
            funcdirs = tmp_funcdirs
        ## reg
        if regdirs:
            regdirs = [ (l,self._psub(x)) for l,x in regdirs.iteritems() ]
            [ 
                self.log.error("registration directory '%s' does not exist" % x)
                for l,x in regdirs if not op.isdir(x)
            ]
        ## nuisance
        if nuidirs:
            tmp_nuidirs = [ (l,glob(self._psub(x))) for l,x in nuidirs.iteritems() ]
            [ 
                self.log.error("nuisance directory '%s' does not exist" % nuidirs[l])
                for l,x in tmp_nuidirs if len(x)==0
            ]
            nuidirs = tmp_nuidirs
        
        ## save
        sinfo['anatdirs'] = anatdirs
        sinfo['segdirs'] = segdirs
        sinfo['regdirs'] = regdirs
        sinfo['funcdirs'] = funcdirs
        sinfo['nuidirs'] = nuidirs
        
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
        
        # seg
        if segdirs:
            self.log.debug("...segmentation")
            for label,segdir in segdirs:
                sinfo[label] = [
                    ('gm', op.join(segdir, "gm.png")), 
                    ('wm', op.join(segdir, "wm.png")), 
                    ('csf', op.join(segdir, "csf.png"))
                ]
                for k,v in sinfo[label]:
                    if not op.isfile(v):
                        self.log.error("segmentation pic '%s' does not exist" % v)
        
        # reg
        if regdirs:
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
        if funcdirs:
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
                    # outliers pic
                    foutlier = op.join(rundir, self._motion_outliers_fname)
                    if not op.isfile(foutlier):
                        self.log.error("Outlier picture '%s' does not exist" % foutlier)
                    tmpinfo['outlier'] = foutlier
                    # save
                    minfo[label].append((run+1,tmpinfo))
        
        # nuisance
        if nuidirs:
            self.log.debug("...nuisance")
            for label,rundirs in nuidirs:
                ninfo[label] = []
                for run,rundir in enumerate(rundirs):
                    tmpinfo = {}
                    tmpinfo['global'] = op.join(rundir, "mask_global.png")
                    tmpinfo['wm'] = op.join(rundir, "mask_wm.png")
                    tmpinfo['csf'] = op.join(rundir, "mask_csf.png")
                    for k,v in tmpinfo.iteritems():
                        if not op.isfile(v):
                            self.log.error("nuisance pic '%s' does not exist" % v)
                    ninfo[label].append((run+1,tmpinfo))
        
        self.log.debug("...saving")
        self.addTemplateContext('subject', None)
        self.subjects.append((subject, {
            'anat': ainfo,
            'seg': sinfo, 
            'motion': minfo, 
            'reg': rinfo, 
            'nuisance': ninfo
        }))
        os.chdir(cwd)
        return
    
    def setData(self, subjects, anatdirs, funcdirs=[], regdirs=[], 
                segdirs=[], nuidirs=[]):
        self.log.info("Setting Data")
        self.subjects = []
        for subject in subjects:
            self.addSubject(subject, anatdirs, funcdirs, regdirs, segdirs, nuidirs)
        return
    

