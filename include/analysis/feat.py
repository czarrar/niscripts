# TODO: split this stuff up

import os, nibabel, re, yaml, tempfile, shutil
import os.path as op
from string import Template
from collections import OrderedDict
from jinja2 import Environment, PackageLoader
from glob import glob
import numpy as np
from copy import deepcopy
from execute import Process
from analysis.base import Base, SubjectBase
from copy import deepcopy
#import nipype.interfaces.fsl as fsl # fsl

# Time
NIFTI_UNITS_SEC=8
NIFTI_UNITS_MSEC=16
NIFTI_UNITS_USEC=24

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class CombineSubject(SubjectBase):
    """
    Combine Different Functional Runs
    
    Examples
    --------
    >>> from feat_subject import CombineSubject
    >>> 
    >>> template_vars = {'base': "/Users/zarrar/Projects/tnetworks/output/${subject}/wm", 
                            'subject': 'tb3417'}
    >>> cf = CombineSubject(2, template_vars)
    >>> 
    >>> infiles = "${base}/run_[1-6]/func_preproc.nii.gz"
    >>> outfunc = "${base}/all_runs/func_preproc.nii.gz"
    >>> cf.setData(infiles, outfunc, mkdir="${base}/all_runs", overwrite=False)
    >>> 
    >>> inconfound = "${base}/motion_all.par"
    >>> outconfound = "${base}/all_runs/confound_evs.1D"
    >>> cf.setDecon(polort="A", censortr="2:37,3:48", inconfound=inconfound, outconfound=outconfound)
    >>> 
    """
    _logname = "combine_subject"
    
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        verbosity : 0=minimal, 1=verbose, 2=debug
        template_vars: template variables for paths
        """
        super(CombineSubject, self).__init__(*args, **kwargs)
        self._isset_data = False
        self._isset_decon = False
        self.log.debug("Starting CombineSubject")
    
    def fromDict(self, config):
        self.check_req(config, ["data", "decon"])
        
        data = config.pop("data")
        self.check_req(data, ["infiles", "outfunc"])
        self.setData(**data)
        
        motion = config.pop("motion", None)
        if motion:
            self.check_req(motion, ["infiles", "outfile"])
            self.setMotion(**motion)
        
        return
    
    def compile(self):
        pass
    
    def run(self):
        if not self._isset_data or not self._isset_decon:
            self.log.fatal("Data or 3dDeconvolve must be setup before running")
        
        # Combine Functionals
        cmd = "fslmerge -t %s %s" % (self.outfunc, " ".join(self.infiles))
        self.log.command(cmd, cwd=op.dirname(self.outfunc))
        
        return
    
    def setData(self, infiles, outfunc, runfile=None, mkdir=None, overwrite=False):
        """Set inputs and outputs of functional data"""
        self.log.debug("Setting Data")
        
        self.infiles = super(CombineSubject, self).getInputs(infiles, runfile, mkdir)
        
        # Set outfunc
        outfunc = self._substitute(outfunc)
        if op.isfile(outfunc):
            if overwrite:
                self.log.warning("removing output '%s'" % outfunc)
                os.remove(outfunc)
            else:
                self.log.error("Output '%s' already exists...not overwriting" % outfunc)
        self.outfunc = outfunc
        
        self._isset_data = True
        
        return
    

class ResDeconSubject(SubjectBase):
    _logname = "res_decon_subject"
    
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        verbosity : 0=minimal, 1=verbose, 2=debug
        template_vars: template variables for paths
        """
        super(ResDeconSubject, self).__init__(*args, **kwargs)
        self._isset_data = False
        self._isset_options = False
        self.log.debug("Starting ResDeconSubject")
    
    def fromDict(self, config):
        self.check_req(config, ["data", "options"])
        
        data = config.pop("data")
        self.check_req(data, ["infiles"])
        self.setData(**data)
        
        motion = config.pop("motion", None)
        if motion:
            self.check_req(motion, ["infiles", "outfile"])
            self.setMotion(**motion)
        
        options = config.pop("options")
        self.setOptions(**options)
        
        return
    
    def compile(self):
        return
    
    def run(self):
        self.compile()
        
        if not self._isset_data or not self._isset_options:
            self.log.critical("Must specify data or options")
        
        # decon
        if self.decon_opts:
            if self.outfunc:
                cwd = op.dirname(self.outfunc)
            elif self.outmat:
                cwd = op.dirname(self.outmat)
            self.log.command(" ".join(self.decon_opts), cwd=cwd)
        else:
            raise Exception("todo")
        
        # confound file
        if self.outconfound:
            x = np.loadtxt(self.outmat)
            np.savetxt(self.outconfound, x, fmt="%.9f")
        if self.outmat_istmp:
            os.remove(self.outmat)
        
        # add back mean
        if self.outfunc:
            self.log.info("Adding mean to residuals")
            
            self.log.debug("Merging original functional data")
            tmpfunc = op.join(op.dirname(self.outfunc), "tmp.nii.gz")
            if op.isfile(tmpfunc):
                self.log.warning("Removing tmpfunc %s" % tmpfunc)
                os.remove(tmpfunc)
            cmd = "fslmerge -t %s %s" % (tmpfunc, " ".join(self.infiles))
            self.log.command(cmd, cwd=op.dirname(self.outfunc))
            
            self.log.debug("Getting mean")
            meanfunc = op.join(op.dirname(self.outfunc), "func_mean.nii.gz")
            if op.isfile(meanfunc):
                self.log.warning("Removing meanfunc %s" % meanfunc)
                os.remove(meanfunc)
            cmd = "fslmaths %s -Tmean %s " % (tmpfunc, meanfunc)
            self.log.command(cmd, cwd=op.dirname(self.outfunc))
            
            self.log.debug("Adding mean")
            tmpfunc = op.join(op.dirname(self.outfunc), "tmp.nii.gz")
            if op.isfile(tmpfunc):
                self.log.warning("Removing tmpfunc %s" % tmpfunc)
                os.remove(tmpfunc)
            cmd = "3dcalc -a %s -b %s -expr 'a+b' -prefix %s" % (self.outfunc, meanfunc, tmpfunc)
            self.log.command(cmd, cwd=op.dirname(self.outfunc))
            self.log.info("Moving %s => %s" % (tmpfunc, self.outfunc))
            shutil.move(tmpfunc, self.outfunc)
        
        return
    
    def setData(self, infiles, outfunc="", outmat="", outconfound="", runfile=None, mkdir=None, 
                overwrite=False):
        """Set inputs and outputs of functional data"""
        self.log.debug("Setting Data")
        
        self.infiles = super(ResDeconSubject, self).getInputs(infiles, runfile, mkdir)
        
        # Set outfunc
        outfunc = self._substitute(outfunc)
        if outfunc and op.isfile(outfunc):
            if overwrite:
                self.log.warning("removing output '%s'" % outfunc)
                os.remove(outfunc)
            else:
                self.log.error("Output '%s' already exists...not overwriting" % outfunc)
        self.outfunc = outfunc
        
        # Set outmat
        outmat = self._substitute(outmat)
        if outmat and op.isfile(outmat):
            if overwrite:
                self.log.warning("removing output '%s'" % outmat)
                os.remove(outmat)
            else:
                self.log.error("Output '%s' already exists...not overwriting" % outmat)
        self.outmat = outmat
        
        # Outconfound
        self.outconfound = self._substitute(outconfound)
        if not outmat and outconfound:
            self.outmat_istmp = True
            self.outmat = tempfile.mktemp(".1D")
        else:
            self.outmat_istmp = False
        
        # check
        if not outfunc and not outmat:
            self.log.critical("Must specify either outfunc or outmat")
        
        self._isset_data = True
        return
    
    def setOptions(self, **kwargs):
        if not self._isset_data:
            self.log.fatal("Must set data before options")
        
        decon_opts = []
        
        # Have Anything?
        if kwargs == {}:
            self.decon_opts = decon_opts
            return
        
        # Defaults
        kwargs.setdefault('input', " ".join(self.infiles))
        
        if self.outmat:
            kwargs.setdefault('x1D', self.outmat)
        
        if self.outfunc:
            kwargs.setdefault('errts', self.outfunc)
            kwargs.setdefault('nocout', True)
            kwargs.setdefault('nobucket', True)
        else:
            kwargs.setdefault('x1D_stop', True)
        
        # Set user kwargs
        decon_opts.append("3dDeconvolve")
        for k,v in kwargs.iteritems():
            if isinstance(v, bool) and v == True:
                decon_opts.append("-%s" % k)
            else:
                decon_opts.append("-%s %s" % (k,v))
        
        self.decon_opts = decon_opts
        
        self._isset_options = True
        return
    

class TsSubject(SubjectBase):
    def __init__(self, *args, **kwargs):
        super(TsSubject, self).__init__(*args, **kwargs)
        self.log.info("Starting TsSubject")
        self._isset_data = False
        self.data_opts = []
        self._isset_options = False
        self.option_opts = []
        self.cmd = ""
        return
    
    def fromDict(self, config):
        self.check_req(config, ["data"])
        
        data = config.pop("data")
        self.check_req(data, ["infunc", "outts"])
        self.setData(**data)
        
        options = config.pop("options", None)
        if options:
            self.setOptions(**options)
        
        return
    
    def compile(self):
        self.log.info("Compiling")
        if not self._isset_data or not self.data_opts:
            self.log.error("Must set data first")
        cmd_opts = ["fslmeants"] + self.data_opts + self.option_opts
        self.cmd = " ".join(cmd_opts)
        return self.cmd
    
    def run(self):
        self.compile()
        self.log.info("Running")
        
        if not self.cmd:
            self.log.fatal("No command")
        
        if self.dry_run:
            self.log.drycommand(self.cmd)
        else:
            self.log.command(self.cmd)
        
        return
    
    def setData(self, infunc, outts, inmask=None, incoord=None, overwrite=False):
        self.log.info("Setting data")
        data_opts = []
        
        # Input 4D
        infunc = self._substitute(infunc)
        if not op.isfile(infunc):
            self.log.error("Input 4D functional '%s' does not exist" % infunc)
        data_opts.append("-i %s" % infunc)
        
        # ROI
        if inmask is not None and incoord is not None:
            self.log.error("Cannot specify both input mask (inmask) and input coords (incoord)")
        elif inmask:
            inmask = self._substitute(inmask)
            if not op.isfile(inmask):
                self.log.error("Input mask (inmask=%s) does not exist" % inmask)
            data_opts.append("-m %s" % inmask)
        elif incoord:
            if not isinstance(incoord, list) or len(incoord) != 3:
                self.log.error("Input coordinates (incoord=%s) must be a list of length 3" % 
                                    incoord)
            data_opts.append("-c %s" % " ".join(incoord))
        else:
            self.log.error("Must specify either input mask (inmask) or input coords (incoord)")
        
        # Output
        outts = self._substitute(outts)
        if op.isfile(outts):
            if overwrite:
                self.log.warning("Removing output '%s'" % outts)
                os.remove(outts)
            else:
                self.log.error("Output '%s' already exists" % outts)
        if not op.isdir(op.dirname(outts)):
            self.log.info("Creating base directory '%s' of output" % op.dirname(outts))
            os.mkdir(op.dirname(outts))
        data_opts.append("-o %s" % outts)
        
        self.data_opts = data_opts
        self._isset_data = True
        return 
    
    def setOptions(self, **kwargs):
        self.log.info("Setting options")
        option_opts = []
        
        choices = {
            "usemm": 'opt', 
            "showall": 'opt', 
            "eig": 'opt', 
            "order": 'optarg', 
            "no_bin": 'opt', 
            "label": 'optfile', 
            "transpose": 'opt', 
            "verbose": 'opt'
        }
        
        for k,v in kwargs.iteritems():
            if k not in choices:
                self.log.error("Bad option '%s'; must be one of %s" % 
                                    (k, ", ".join(choices.keys())))          
            if choices[k] == "opt":
                if not isinstance(v, bool):
                    self.log.error("Option %s must have a value of True/False" % v)
                if v:
                    option_opts.append("--%s" % k)
            elif choices[k] == "optarg":
                option_opts.append("--%s %s" % (k,v))
            elif choices[k] == "optfile":
                if not op.isfile(v):
                    self.log.error("Couldn't find file '%s' for option '%s'" % (v,k))
                option_opts.append("--%s %s" % (k,v))
            else:
                raise Exception("todo")
        
        self.option_opts = option_opts
        self._isset_options = True
        return
    

class FsfSubject(SubjectBase):
    """
    Examples
    --------
    
    >>> from feat_subject import FsfSubject
    >>> fsf = FsfSubject(verbosity=2)
    >>> 
    >>> base = "/Users/zarrar/Projects/tnetworks/output/tb3417/wm/run_1/"
    >>> 
    >>> infile = base + "func_preproc.nii.gz"
    >>> outdir = base + "ztest.feat"
    >>> fsf.setData(infile=infile, outdir=outdir, highpass_filter=128)
    >>> 
    >>> confoundev_file = base + "motion.par"
    >>> fsf.setStats(confoundev_file=confoundev_file)
    >>> 
    >>> fsf.setPostStats()
    >>> 
    >>> evfname = base + "ev1_test.txt"
    >>> fsf.addEV("mem", evfname, "double-gamma", tempderiv=True)
    >>> evfname = base + "ev2_test.txt"
    >>> fsf.addEV("del1", evfname, "double-gamma", tempderiv=True)
    >>> 
    >>> fsf.addContrast("memoranda", "+1*mem")
    >>> fsf.addContrast("delay1", "+1*del1")
    >>> fsf.addContrast("memorand>delay1", "+0.5*mem -0.5*del1")
    >>> fsf.addContrast("memorand+delay1", "+0.5*mem +0.5*del1")
    >>> 
    >>> print fsf.compile()
    >>> fsf.write()
    >>> 
    
    """
    
    _logname = "fsf_subject"
    
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        verbosity : 0=minimal, 1=verbose, 2=debug
        template_vars: template variables for paths
        """
        super(FsfSubject, self).__init__(*args, **kwargs)
        
        self.env = Environment(loader=PackageLoader('analysis', 'templates'))
        
        self.bytrial_evnames = []
        self.fsf_context = {}
        self.evs = OrderedDict()
        self.contrasts = OrderedDict()
        self.real_contrasts = OrderedDict()
        self.ftests = []
        
        isset_names = ['data', 'stats', 'poststats']
        self._isset = dict.fromkeys(isset_names, False)
        self._is_compiled = False
        
        self.fsf = []
        
        self.log.debug("Starting FsfSubject")
        return
    
    def fromDict(self, config):
        self.check_req(config, ["data", "stats", "poststats"])
        
        # Data
        data = config.pop("data")
        self.check_req(data, ["infile", "outdir", "outfsf", "highpass_filter"])
        self.setData(**data)
        
        # Stats
        stats = config.pop("stats")
        self.setStats(**stats)
        
        # Post-Stats
        poststats = config.pop("poststats")
        if poststats:
            self.setPostStats(**poststats)
        else:
            self.setPostStats(do_poststats=False)
        
        if len(config) > 0:
            self.log.error("Stuff left-over from config: %s" % config)
        
        return
    
    def run(self, check=True):
        fsf = self.compile()
        self.write(check=check)
        return fsf
    
    def write(self, check=True):
        sfsf = self._remove_extra_stuff("\n".join(self.fsf))
        f = open(self.outfsf, 'w')
        f.write(sfsf)
        f.close()
        if check:
            if 'confoundev_file' in self.fsf_context:
                conf = " " + self.fsf_context['confoundev_file']
            else:
                conf = ""
            self.log.command("feat_model %s%s" % (op.splitext(self.outfsf)[0], conf))
        return
    
    def compile(self, redo=False):
        self.log.debug("Compiling fsf...")
        
        # Checks
        self.log.debug("checks")
        if self._is_compiled and not redo:
            raise Exception("Can't compile FEAT model more than once")
        for name,isset in self._isset.iteritems():
            if not isset:
                self.log.error("Have not set %s" % name)
        if len(self.evs) == 0:
            self.log.error("You must specify at least one EV")
        if len(self.contrasts) == 0:
            self.log.error("You must specify at least one contrast")
        
        # Get 'real' contrasts
        self.log.debug("get real contrasts")
        if not redo:
            self._setRealContrasts()
        
        # Double-Check
        self.log.debug("double checks")
        ## contrasts
        con_nref = len(self.contrasts[self.contrasts.keys()[0]])
        if con_nref != len(self.evs):
            raise Exception("todo")
        for k,v in self.contrasts.iteritems():
            if len(v) != con_nref:
                self.log.error("Inconsistent number of contrasts (%i) for %s" % (len(v), k))
        ## real contrasts
        con_real_nref = len(self.real_contrasts[self.real_contrasts.keys()[0]])
        for k,v in self.real_contrasts.iteritems():
            if len(v) != con_real_nref:
                self.log.error("Inconsistent number of real contrasts (%i) for %s" % (len(v), k))        
        ## f-tests
        if self.ftests:
            nref = len(self.ftests[0])
            for v in self.ftests:
                if len(v) != nref:
                    self.log.error("Inconsistent number of ftests (%i vs %i)" % (len(v), nref))
            
        # Set nums
        self.log.debug("set number of evs, cons, and ftests")
        nevs = con_nref
        nevs_real = con_real_nref
        ncons = len(self.contrasts)
        nftests = len(self.ftests)
        self.fsf_context.update({
            'num_evs': nevs,
            'num_evs_real': nevs_real,
            'num_cons': ncons,
            'num_fcons': nftests
        })
        
        # Gather fsf        
        fsf = []
        ## header
        self.log.debug("get header")
        fsf_part = self._render('subject_header', **self.fsf_context)
        fsf.append(fsf_part)
        ## evs
        self.log.debug("get evs")
        fsf_part = self._render('subject_evs', evs=self.evs)
        fsf.append(fsf_part)
        ## cons
        self.log.debug("get cons")
        fsf_part = self._render('subject_contrasts', con_type="real", 
                                contrasts=self.real_contrasts, ftests=self.ftests, 
                                nevs=nevs_real, ncons=ncons, nftests=nftests)
        fsf.append(fsf_part)
        fsf_part = self._render('subject_contrasts', con_type="orig", contrasts=self.contrasts, 
                                ftests=self.ftests, nevs=nevs, ncons=ncons, nftests=nftests)
        fsf.append(fsf_part)
        ## con masks
        self.log.debug("get con masking")
        fsf_part = self._render('subject_contrast_masking', n=len(self.contrasts)+len(self.ftests))
        fsf.append(fsf_part)
        ## end
        self.log.debug("get non-gui end")
        fsf_part = self._render('subject_nongui', **self.fsf_context)
        fsf.append(fsf_part)
        
        self._is_compiled = True
        self.fsf = fsf
        
        return self._remove_extra_stuff("\n".join(fsf))
    
    def recompile(self):
        return self.compile(redo=True)
    
    def _remove_extra_stuff(self, a):
        b = re.sub("\n\ +", "\n", a)
        c = re.sub("\n\n+", "\n\n", b)
        return c
    
    def _render(self, fname, **context):
        template = self.env.get_template(fname + '.jinja')
        result = template.render(**context)
        return self._remove_extra_stuff(result)
    
    def setData(self, infile, outdir, outfsf, highpass_filter, 
                hp_yn=False, tr=None, overwrite=0):
        infile = self._substitute(infile)
        outdir = self._substitute(outdir)
        outfsf = self._substitute(outfsf)
        
        # infile
        gpaths = glob(infile)
        if len(gpaths) == 1:
            infile = gpaths[0]
        elif len(gpaths) > 1:
            self.log.error("Cannot have more than one infile %s from %s" % (gpaths, infile))
        else:
            self.log.error("Did not find infile %s" % infile)
        if not op.isfile(infile):
            self.log.error("Could not find input file %s" % infile)
        
        self.fsf_context['func_file'] = infile
        self.fsf_context['outputdir'] = outdir
        self.fsf_context['highpass_yn'] = hp_yn
        self.fsf_context['highpass_filter'] = highpass_filter
        
        func = nibabel.load(infile)
        if len(func.shape) != 4:
            self.log.error("Input file '%s' should have 4 dimensions" % infile)
        self.fsf_context['num_vols'] = func.shape[3]
        
        if not tr:
            hdr = func.get_header()
            tr = hdr['pixdim'][4]
            if hdr['xyzt_units'] > NIFTI_UNITS_USEC:
                tr /= 1e6
            elif hdr['xyzt_units'] > NIFTI_UNITS_MSEC:
                tr /= 1e3
        self.fsf_context['tr'] = tr
        
        if overwrite == 0 and op.isdir(outdir):
            self.log.error("Output directory '%s' already exists" % outdir)
        elif overwrite > 0:
            self.fsf_context['overwrite'] = overwrite - 1
        
        if op.isdir(outdir):
            if overwrite:
                self.log.warning("Output directory already exists, will overwrite")
            else:
                self.log.warning("Output directory already exists, instead of overwriting will create a new output directory")
        dir_outdir = op.dirname(outdir)
        if not op.isdir(dir_outdir):
            self.log.info("Creating directory '%s' for feat output" % dir_outdir)
            os.mkdir(dir_outdir)
        
        dir_outfsf = op.dirname(outfsf)
        if not op.isdir(dir_outfsf):
            self.log.info("Creating directory '%s' for fsf output" % dir_outfsf)
            os.mkdir(dir_outfsf)
        self.outfsf = op.splitext(outfsf)[0] + ".fsf"
        
        self._isset['data'] = True
        return
    
    def setStats(self, prewhiten=True, confoundev_file="", evs=[], contrasts=[]):
        self.fsf_context['prewhiten'] = prewhiten
        
        if confoundev_file:
            confoundev_file = self._substitute(confoundev_file)
            if not op.isfile(confoundev_file):
                self.log.error("Could not find confound EV file %s" % confoundev_file)
            self.fsf_context['confoundev_yn'] = True
        else:
            self.fsf_context['confoundev_yn'] = False
        self.fsf_context['confoundev_file'] = confoundev_file
        
        for ev in evs:
            ev_name,ev_opts = ev.items()[0]
            self.addEV(ev_name, **ev_opts)
        for contrast in contrasts:
            con_name,con_val = contrast.items()[0]
            self.addContrast(con_name, con_val)
        
        self._isset['stats'] = True
        return
    
    def setPostStats(self, thresh_type="cluster", vox_thresh=2.3, clust_thresh=0.05, prethreshold_mask="", do_poststats=True):
        self.poststats = do_poststats
        
        thresh_type_choices = {
            'none':         0,
            'uncorrected':  1,
            'voxel':        2,
            'cluster':      3
        }
        try:
            self.fsf_context['thresh_type'] = thresh_type_choices[thresh_type]
        except KeyError:
            self.log.error("Unrecognized threshold type '%s'" % thresh_type)
        
        prethreshold_mask = self._substitute(prethreshold_mask)
        self.fsf_context['prethreshold_mask'] = prethreshold_mask
        if prethreshold_mask and not op.isfile(prethreshold_mask):
                self.log.error("Could not find pre-threshold mask %s" % prethreshold_mask)
        
        self._isset['poststats'] = True
        return
    
    
    def addEV(self, name, fname, convolution, tempfilt=False, tempderiv=False, 
              bytrial=None, bycolumn=None, **opts):
        """
        Adds an Expanatory Variable (EV) to the evs list.
        See FEAT help for more details on inputs.
        
        Parameters
        ----------
        name : name of EV
        fname : timing file with 1 or 3 columns
        convolution: 'none', 'gamma', or 'double-gamma' are the options
        tempfilt: do temporal filtering?
        tempderiv: get temporal derivative?
        **opts: any additional options
        
        """
        
        waveform_choices = {
            'custom1': 2,
            'custom3': 3,
            'empty': 10     # NOT used
        }
        convolution_choices = {
            'none': 0,
            'gamma': 2,
            'double-gamma': 3,
            'basis': 7      # NOT used
        }
        
        fname = self._substitute(fname)
        
        # Checks
        if self.contrasts:
            raise Exception("Cannot add EVs after contrasts has been set")
        if name in self.evs:
            self.log.error("Can't add EV '%s', since it already exists" % name)
        if not op.isfile(fname):
            self.log.error("File '%s' for EV '%s' does not exist" % (fname, name))
        if bytrial is not None and bycolumn is not None:
            self.log.error("Cannot specify both bytrial and bycolumn")
        
        # Waveform Shape (1 column or 3 column file)
        ncols = self._file_ncols(fname)
        if bycolumn is not None:
            waveform = waveform_choices['custom1']
        else:
            if ncols == 1:
                waveform = waveform_choices['custom1']
            elif ncols == 3:
                waveform = waveform_choices['custom3']
            if ncols != 1 and ncols != 3:
                self.log.error("File %s for EV %s must have 1 or 3 columns, instead %i were found" % (fname, name, ncols))
        
        # Convolution
        try:
            orig_conv = convolution
            convolution = convolution_choices[convolution]
        except KeyError as k:
            self.log.error("Convolution '%s' not recognized" % convolution)
        
        # Calculate for individual trials?
        if bytrial is not None:
            self.log.info("ev %s will be split into individual trials" % name)
            if not isinstance(bytrial, list):
                self.log.fatal("bytrial must be a list of 2 arguments [WHICH_TRIALS, TMP_DIRECTORY]")
            if ncols != 3:
                self.log.fatal("if calculating betas for each individual trial, then input file must be 3" \
                               " column format")
            f = open(fname, 'r')
            lines = f.readlines()
            lines = [ l.strip() for l in lines if l.strip() ]
            ntrials = len(lines)
            f.close()
            which_trials, basedir = tuple(bytrial)
            which_trials = self._substitute(which_trials)
            basedir = self._substitute(basedir)
            if not op.isdir(basedir):
                self.log.info("Creating directory %s" % basedir)
                os.mkdir(basedir)
            
            if which_trials == 'all':
                for i in xrange(ntrials):
                    sname = "%s_trial%04i" % (name, i+1)
                    fn = op.join(basedir, "%s.txt" % sname)
                    # create new ev files
                    f = open(fn, 'w')
                    f.write(lines[i] + "\n")
                    f.close()
                    # call ev function
                    self.addEV(sname, fn, orig_conv, tempfilt, tempderiv, 
                               bytrial=None, **opts)
            else:
                i = int(which_trials)
                sname = "%s_trial_%04i" % (name, i)
                fn = op.join(basedir, "%s.txt" % sname)
                f = open(fn, 'w')
                f.write(lines[i-1] + "\n")
                f.close()
                self.addEV(sname, fn, orig_conv, tempfilt, tempderiv, 
                           bytrial=None, **opts)
            return
        elif bycolumn is not None:
            self.log.info("ev %s will be split into individual columns" % name)
            if not isinstance(bycolumn, list):
                self.log.fatal("bycolumn must be a list of 2 arguments [WHICH_COLS, TMP_DIRECTORY]")
            x = np.loadtxt(fname)
            which_cols, basedir = tuple(bycolumn)
            which_cols = self._substitute(which_cols)
            basedir = self._substitute(basedir)
            if not op.isdir(basedir):
                self.log.info("Creating directory %s" % basedir)
                os.mkdir(basedir)
            
            if which_cols == 'all':
                for i in xrange(ncols):
                    sname = "%s_col%04i" % (name, i+1)
                    fn = op.join(basedir, "%s.txt" % sname)
                    np.savetxt(fn, x[:,i], fmt="%f")
                    self.addEV(sname, fn, orig_conv, tempfilt, tempderiv, 
                               bycolumn=None, **opts)
            else:
                i = int(which_cols)
                sname = "%s_col%04i" % (name, i)
                fn = op.join(basedir, "%s.txt" % sname)
                np.savetxt(fn, x[:,i-1], fmt="%f")
                self.addEV(sname, fn, orig_conv, tempfilt, tempderiv, 
                           bycolumn=None, **opts)
            return
        else:
            ev = {
                'waveform': waveform,
                'convolution': convolution,
                'fname': fname,
                'tempfilt': tempfilt,
                'tempderiv': tempderiv
            }
            ev.update(opts)
            self.evs[name] = ev
        
            self._isset['ev'] = True
            return
    
    def addContrast(self, name, con):
        """
        Add contrast to list of contrasts.
        
        Parameters
        ----------
        con : can be a string or list...
        
        """
        
        if name in self.contrasts:
            self.log.error("Contrast %s already exists" % name)
        if self.ftests:
            self.log.fatal("Cannot add contrasts after setting F-tests")
        if not self.evs:
            self.log.fatal("Add all your EVs before adding any contrasts")
        if self.real_contrasts:
            self.log.fatal("Cannot add contrasts after setting real contrasts")
        
        if isinstance(con, list) and len(con) != len(self.evs):
            self.log.error("Contrast %s list must have a len=%i" % (name, len(self.evs)))
        elif isinstance(con, str):
            final_con = [ 0 for x in xrange(len(self.evs)) ]
            elems = re.split("[ \t]+", con)
            ev_names = self.evs.keys()
            for e in elems:
                val, ev_name = tuple(e.split("*"))
                try:
                    i = ev_names.index(ev_name)
                    final_con[i] = float(val)
                except ValueError:
                    self.log.error("Could not find EV %s for contrast %s" % (ev_name, name))
            con = final_con
        
        self.contrasts[name] = con
        
        self._isset['contrast'] = True
        return
    
    def addFtest(self, ftest):
        if not self.contrasts:
            self.log.fatal("Add all your contrasts before adding any F-tests")
        
        if isinstance(ftest, list) and len(ftest) != len(self.contrasts):
            self.log.error("Ftest must have a len=%i" % len(self.contrasts))
        elif isinstance(ftest, str):
            final_ftest = [ 0 for x in xrange(len(self.contrasts)) ]
            elems = re.split("[, \t]+", ftest)
            con_names = self.contrasts.keys()
            for e in elems:
                try:
                    i = con_names.index(e)
                    final_ftest[i] = 1
                except ValueError:
                    self.log.error("Could not find contrast %s for given ftest" % e)
            ftest = final_ftest
        
        self.ftests.append(ftest)
        
        self._isset['ftest'] = True
        return
    
    def _setRealContrasts(self):
        if self.real_contrasts:
            self.log.fatal("Real contrasts has already been set")
        if not self.contrasts:
            self.log.fatal("Orig contrasts have not yet been created")
        
        ev_names = self.evs.keys()
        for i,con in enumerate(self.contrasts):
            con_vals = self.contrasts[con]
            real_cons = []
            for j,v in enumerate(con_vals):
                real_cons.append(v)
                k = ev_names[j]
                if self.evs[k]['tempderiv']:
                    real_cons.append(0)
            self.real_contrasts[con] = real_cons
        
        return
    

class BetaSeriesSubject(SubjectBase):
    
    _logname = "beta_series_subject"
    
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        verbosity : 0=minimal, 1=verbose, 2=debug
        template_vars: template variables for paths
        """
        print 'hey'
        self.verbosity = args[0]
        super(BetaSeriesSubject, self).__init__(*args, **kwargs)
        self.log.debug("Starting BetaSeriesSubject")
        return
    
    def fromDict(self, config):
        self.check_req(config, ["data", "stats", "bs"])
        # Beta-Series
        bs = config.pop("bs")
        self.setBetaSeries(**bs)
        # Save the rest of the config
        self.setConfig(config)
        return
    
    def compile(self):
        return
    
    def run(self):
        if self._isrun:
            self.log.fatal("beta-series already run")
        self.compile()
        
        # TODO: automask? with np.std along the cols of y?
        self.log.info("loading mask")
        mask = nibabel.load(self.maskfile)
        m = mask.get_data()
        m = m.reshape(np.product(m.shape))
        inds = m.nonzero()[0]
        
        self.log.info("loading data")
        func = nibabel.load(self.infile)
        y = func.get_data()
        # reshape to 2D matrix with rows as time-points and cols as voxels
        nvoxs = np.product(y.shape[:-1])
        ntpts = y.shape[-1]
        y = y.reshape((nvoxs, ntpts)).transpose()
        
        self.log.info("masking data")
        y = y[:,inds]
        
        self.log.info("creating empty beta-series image")
        bseries = np.zeros((nvoxs, self.ntrials), dtype=np.float32)
        
        self.log.info("work time")
        for i in self.trials:
            ii = i - 1
            self.log.title("Trial #%i", i)
            
            # fsf
            self.log.subtitle("creating design matrix")
            tmp_context = deepcopy(self.template_context)
            tmp_context['data']['outdir'] = self.fsf_outdir
            outprefix = "%s_trial%04i" % (self.fsf_prefix, i)
            tmp_context['data']['outfsf'] = outprefix + ".fsf"
            tmp_context['bs_trial'] = i
            fsf = FsfSubject(self.verbosity, tmp_context, 
                             dry_run=self.dry_run, logger=self.log)
            fsf.fromDict(deepcopy(self.config))
            fsf.run()
            
            # regression
            self.log.subtitle("running regression")
            ## get the generated mat file
            outmat = outprefix + ".mat"
            X = np.loadtxt(outtxt, skiprows=5)
            ## calculate betas (code borrowed from statmodels)
            pinv_X = np.linalg.pinv(X)
            betas = np.dot(pinv_X, y)
            ## save
            bseries[:,inds[ii]] = betas[self.ev_index,:]
        
        # save as nifti
        self.log.info("saving beta-series")
        bseries.shape = tuple(list(func.shape[:3]) + [self.ntrials])
        outfname = op.join(self.outdir, "betaseries.nii.gz")
        hdr = func.get_header()
        hdr.set_data_shape(bseries.shape)
        nii_bseries = nibabel.Nifti1Image(bseries, func.get_affine(), hdr)
        nii_bseries.to_filename(outfname)
        
        # save name in text file
        f = file(op.join(self.outdir, "name.txt"), 'w')
        f.write(self.name, "\n")
        f.close()
        
        return
    
    def setBetaSeries(self, name, fname, trials, outdir, maskfile, 
                      outtype='coefs', overwrite=False):
        """Set the EV that will be used to create the beta-series
        
        Parameters
        ----------
        name : name of the EV
        fname : filename with onsets of EV (must be 3 column format)
        trials : can be 'all' for all trials, specified like [1,2,3,4], range like '1:3'
        """
        
        self.log.debug("setting beta-series options")
        
        # Mask File
        self.maskfile = self._substitute(maskfile)
        if not op.isfile(self.maskfile):
            self.log.error("couldn't find maskfile '%s'", self.maskfile)
        
        # Name
        self.name = name
        
        # Total Number of Trials
        fname = self._substitute(fname)
        f = open(fname, 'r')
        lines = f.readlines()
        lines = [ l.strip() for l in lines if l.strip() ]
        ntrials = len(lines)
        f.close()
        
        # Trials
        if isinstance(trials, str):
            l = trials.find(":")
            start = 1
            end = ntrials
            step = 1
            if trials == 'all':
                pass
            elif trials == 'even':
                start = 2
                step = 2
            elif trials == 'odd':
                step = 2
            elif l == 0:
                start = 1
                end = int(trials[1:])
            elif l == (len(tmp)-1):
                start = int(tmp[:l])
                end = ntrials
            elif l != -1:
                start = int(tmp[:l])
                end = int(tmp[(l+1):])
            else:
                self.log.error("cannot parse beta-series trials '%s'" % trials)
            trials = range(start, end+1, step)
        elif isinstance(trials, list):
            trials = [ int(x) for x in trials ]
        else:
            self.log.error("cannot parse beta-series trial '%s'" % trials)
        
        if len(trials) < 2:
            self.log.error("must have at least 2 trials for beta-series")
        
        self.trials = trials
        self.ntrials = len(self.trials)
        
        # Output
        outdir = self._substitute(outdir)
        if op.isdir(outdir):
            if overwrite:
                self.log.warning("Removing output directory '%s'" % outdir)
                shutil.rmtree(outdir)
            else:
                self.log.warning("Output directory '%s' already exists" % outdir)
        else:
            self.log.debug("creating beta-series output directory '%s'" % 
                           outdir)
            os.mkdir(outdir)
        
        dir_outdir = op.dirname(outdir)
        if not op.isdir(dir_outdir):
            self.log.debug("creating directory '%s' for beta-series output" % 
                            dir_outdir)
            os.mkdir(dir_outdir)
        
        self.outdir = outdir
        
        # Output Type
        outtype_choices = ['coefs', 'tstats']
        if outtype not in outtype_choices:
            self.log.error("unrecognized outtype '%s'" % outtype)
        self.outtype = outtype
        
        # Stuff for fsf
        self.fsf_outdir = op.join(self.outdir, "dont_run.feat")
        fsfdir = op.join(self.outdir, "fsfs")
        if not op.isdir(fsfdir):
            self.log.debug("creating main fsf dir")
            os.mkdir(fsfdir)
        self.fsf_prefix = op.join(fsfdir, "design")
        
        self._isset_bs = True
        return
    
    def setConfig(config):
        self.log.debug("config in beta-series")
        if not self._isset_bs:
            self.log.fatal("must set beta-series first")
        self.check_req(config, ["data", "stats"])
        if self.name not in config['stats']['evs']:
            self.log.error("couldn't find beta-series '%s' in EVs", self.name)
        self.ev_index = config['stats']['evs'].keys().index(self.name)
        infile = self._substitute(config["data"]["infile"])
        if not op.isfile(infile):
            self.log.error("input file '%s' does not exist", infile)
        self.infile = infile
        self.config = config
        return
    

class FsfInfo(object):
    def __init__(self, fname):
        """docstring for __init__"""
        f = file(fname, 'r')
        self.fsf = [ l.strip() for l in f if l.strip() != '' and l[0] != '#' ]
        return
    
    def getOutput(self):
        res = [ re.search("set\ fmri\(.*?\)\ (.*)", l) for l in self.fsf if l.find("outputdir") != -1 ]
        res = [ x.groups()[0].strip('"') for x in res ]
        return res

class FeatSubject(SubjectBase):
    """docstring for FeatSubject"""
    
    _logname = "feat_subject"
    
    def __init__(self, *args, **kwargs):
        super(FeatSubject, self).__init__(*args, **kwargs)
    
    def fromDict(self, config_dict):
        self.check_req(config_dict, ["infsf"])
        self.setData(**config_dict)
    
    def run(self):
        # Run feat
        self.log.info("Running feat")
        cmd = "feat %s" % self.infsf
        self.log.command(cmd)
        
        # Get outputdir
        self.log.debug("Getting output directory")
        fsfinfo = FsfInfo(self.infsf)
        outputdir = fsfinfo.getOutput()[0]
        
        # Soft link regdir
        self.log.info("Soft linking registration directory")
        os.symlink(self.regdir, op.join(outputdir, "reg"))
        
        return
    
    def setData(self, infsf, regdir=None):
        infsf = self._substitute(infsf)
        if not op.isfile(infsf):
            self.log.error("Can't find input fsf file '%s'" % infsf)
        self.infsf = infsf
        
        if regdir:
            regdir = self._substitute(regdir)
            if not op.isdir(regdir):
                self.log.error("Can't find input registration directory '%s'" % infsf)
            self.regdir = regdir
        
        return

class ApplyRegSubject(SubjectBase):
    """Apply standard registration to image in functional space
    """
    
    _logname = "apply_reg_subject"
    
    def __init__(self, *args, **kwargs):
        super(ApplyRegSubject, self).__init__(*args, **kwargs)
        self.log.info("Starting ApplyRegSubject")
        self.infiles = None
        self.regdir = None
        self.outfiles = None
        self.regtype = None
        self.many_cmd_opts = None
        return
    
    def fromDict(self, config):
        self.check_req(config, ["data"])
        
        data = config.pop("data")
        self.setData(**data)
        
        options = config.pop("options", {})
        self.setOptions(**options)
        
        reg_files = config.pop("reg_files", {})
        self.setRegFiles(**reg_files)
        
        for k in config:
            self.log.fatal("Unrecognized option '%s'" % k)
        
        return
    
    def compile(self):
        self.log.info("Compiling")
        if not self._isset_data or not self._isset_reg_files or not self._isset_options:
            self.log.critical("Must set data, registration files, and options before compiling")
        many_cmd_opts = []
        
        if self.regtype == "flirt":
            for infile,outfile in zip(self.infiles, self.outfiles):
                cmd_opts = ["flirt"]
                cmd_opts.append("-in %s" % infile)
                cmd_opts.append("-ref %s" % self.stdfile)
                cmd_opts.append("-applyxfm")
                cmd_opts.append("-init %s" % self.matfile)
                cmd_opts.append("-out %s" % outfile)
                for k,v in self.cmd_kwargs.iteritems():
                    if isinstance(v, bool):
                        if v is True:
                            cmd_opts.append("-%s" % k)
                    else:
                        cmd_opts.append("-%s %s" % (k,v))
                many_cmd_opts.append(cmd_opts)
        elif self.regtype == "fnirt":
            for infile,outfile in zip(self.infiles, self.outfiles):
                cmd_opts = ["applywarp"]
                cmd_opts.append("-i %s" % infile)
                cmd_opts.append("-r %s" % self.stdfile)
                cmd_opts.append("--premat=%s" % self.prematfile)
                cmd_opts.append("-w %s" % self.warpfile)
                cmd_opts.append("-o %s" % outfile)
                for k,v in self.cmd_kwargs.iteritems():
                    if isinstance(v, bool):
                        if v is True:
                            cmd_opts.append("--%s" % k)
                    else:
                        cmd_opts.append("--%s=%s" % (k,v))
                many_cmd_opts.append(cmd_opts)
        
        self.many_cmd_opts = many_cmd_opts
        return
    
    def run(self):
        self.compile()
        self.log.info("Running")
        if not self.many_cmd_opts:
            self.log.fatal("Bad")
        for cmd_opts in self.many_cmd_opts:
            cmd = " ".join(cmd_opts)
            if self.dry_run:
                self.log.drycommand(cmd)
            else:
                self.log.command(cmd)
        return
    
    def setData(self, infiles, regdir, outdir=None, outfiles=None, overwrite=False):
        self.log.info("Setting data")
        
        # Set inputs
        infiles = self._getInputsWorker(infiles)
        self.infiles = infiles
        
        # Set outputs
        if outdir is not None and outfiles is not None:
            self.log.error("Cannot specify both output directory (outdir) and output files" \
                            " (outfiles)")
        if outdir:
            outdir = self._substitute(outdir)
            if not op.isdir(outdir):
                self.log.info("Creating output directory")
                os.mkdir(outdir)
            outfiles = [ op.join(outdir, op.basename(x)) for x in infiles ]
            for infile,outfile in zip(infiles, outfiles):
                if op.dirname(infile) == op.dirname(outfile):
                    self.log.error("Output directory '%s' has same base path as input" % outdir)
        elif outfiles:
            if isinstance(outfiles, str):
                outfiles = [outfiles]
            outfiles = [ self._substitute(x) for x in outfiles ]
            if len(outfiles) != len(infiles):
                self.log.error("# of outputs (%i) not same as # of inputs (%i)" % 
                                    (len(outfiles), len(infiles)))
        else:
            self.log.error("Must specify either output directory (outdir) or output files" \
                            " (outfiles)")
        for outfile in outfiles:
            if op.isfile(outfile):
                if overwrite:
                    self.log.warning("Removing output '%s'" % outfile)
                    os.remove(outfile)
                else:
                    self.log.error("Output '%s' already exists" % outfile)
        self.outfiles = outfiles
        
        # Set regdir
        regdir = self._substitute(regdir)
        self.regdir = regdir
        
        self._isset_data = True
        return
    
    def setRegFiles(self, standard=None, warp=None, premat=None, mat=None):
        """Files found in the reg directory"""
        self.log.info("Setting registration files")
        
        # Can I run this function?
        if not self._isset_data or not self.regdir:
            self.log.error("Must set data before setting registration files")
        if not self._isset_options or not self.regtype:
            self.log.error("Must set options before setting registration files")
        
        # Paths
        if standard is None:
            standard = op.join(self.regdir, "run_2_func2std.nii.gz")
        self.stdfile = self._substitute(standard)
        if warp is None:
            warp = op.join(self.regdir, "highres2standard_warp.nii.gz")
        self.warpfile = self._substitute(warp)
        if premat is None:
            premat = op.join(self.regdir, "example_func2highres.mat")
        self.prematfile = self._substitute(premat)
        if mat is None:
            mat = op.join(self.regdir, "example_func2standard.mat")
        self.matfile = self._substitute(mat)
        
        # Flirt or Fnirt?
        if self.regtype == "auto":
            if op.isfile(self.warpfile):
                self.regtype == "fnirt"
            elif op.isfile(self.matfile):
                self.regtype == "flirt"
            else:
                self.log.error("Couldn't determine registration type of %s" % self.regdir)
        
        # Checks
        if not op.isfile(self.stdfile):
            self.log.error("Standard image '%s' does not exist" % self.stdfile)
        if self.regtype == "flirt":
            if not op.isfile(self.matfile):
                self.log.error("Matrix file '%s' does not exist" % self.matfile)
        if self.regtype == "fnirt":
            if not op.isfile(self.warpfile):
                self.log.error("Warp file '%s' does not exist" % self.warpfile)
            if not op.isfile(self.matfile):
                self.log.error("Matrix file '%s' does not exist" % self.prematfile)
        
        self._isset_reg_files = True
        return
    
    def setOptions(self, regtype="auto", **kwargs):
        self.log.info("Setting options")
        
        regtype_options = ["auto", "flirt", "fnirt"]
        if regtype not in regtype_options:
            self.log.error("Registration type '%s' not recognized." % regtype)
        self.regtype = regtype
        
        self.cmd_kwargs = kwargs
        
        self._isset_options = True
        return
    
        
        
class RegressSubject(SubjectBase):
    
    _logname = "regress_subject"
    
    def __init__(self, *args, **kwargs):
        super(RegressSubject, self).__init__(*args, **kwargs)
        self.log.info("Starting RegressSubject")
        self.cmd = None
    
    def fromDict(self, config_dict):
        self.setVars(config_dict)
        return
    
    def setVars(self, config_dict):
        # Checks
        self.check_req(config_dict, ["in", "out", "design"])
        if "filter" not in config_dict and "filter_all" not in config_dict:
            self.log.critical("Must give filter or filter_all option for RegressSubject")
        elif "filter" in config_dict and "filter_all" in config_dict:
            self.log.critical("Cannot give both filter or filter_all options for RegressSubject")
        
        # Substitute paths
        overwrite = config_dict.pop("overwrite", False)
        choices = ["in", "out", "design", "mask"]
        for k in choices:
            if k in config_dict:
                config_dict[k] = self._substitute(config_dict[k])
        
        # Check paths
        if not op.isfile(config_dict["in"]):
            self.log.error("Cannot find --in %s" % config_dict["in"])
        if not op.isfile(config_dict["design"]):
            self.log.error("Cannot find --design %s" % config_dict["design"])
        if not op.isfile(config_dict["mask"]):
            self.log.error("Cannot find --mask %s" % config_dict["mask"])
        if op.isfile(config_dict["out"]):
            if overwrite:
                self.log.warning("Removing previous output %s" % config_dict["out"])
                os.remove(config_dict["out"])
            else:
                self.log.error("Output %s already exists" % config_dict["out"])
        
        self.cmd_opts = config_dict
        return
    
    def compile(self):
        self.log.info("Compiling")
        #filt = fsl.FilterRegressor(**self.cmd_opts)
        #cmd = filt.cmdline
        if "filter_all" in self.cmd_opts:
            p = Process("grep /NumWaves %s" % self.cmd_opts["design"]) | Process("awk '{print $2}'")
            if p.retcode != 0:
                self.log.error("Could not get the number of columns in %s" % 
                                    self.cmd_opts["design"])
            ncols = int(p.stdout)
            del self.cmd_opts["filter_all"]
            self.cmd_opts["filter"] = ",".join([ str(x+1) for x in xrange(ncols) ])
        
        cmd = ["fsl_regfilt"]
        for k,v in self.cmd_opts.iteritems():
            if len(k) == 1:
                pre = "-"
            else:
                pre = "--"
            if isinstance(v, bool) and v is True:
                cmd.append("%s%s" % (pre, v))
            elif len(k) == 1:
                cmd.append("%s%s %s" % (pre, k, v))
            else:
                cmd.append("%s%s=%s" % (pre, k, v))
        
        self.cmd = " ".join(cmd)
        return self.cmd
    
    def run(self):
        self.compile()
        self.log.info("Running fsl_regfilt")
        self.log.command(self.cmd, cwd=op.dirname(self.cmd_opts["out"]))
    

class FsfGroup(Base):
    _logname = "fsf_group"
    
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        verbosity : 0=minimal, 1=verbose, 2=debug
        template_vars: template variables for paths
        """
        super(FsfGroup, self).__init__(*args, **kwargs)
        
        self.env = Environment(loader=PackageLoader('analysis', 'templates'))
                
        self.fsf_context = {}
        self.inputs = []
        self.groups = []
        self.evs = OrderedDict()
        self.contrasts = OrderedDict()
        self.ftests = []
        
        isset_names = ['data', 'evs', 'contrasts', 'stats', 'poststats']
        self._isset = dict.fromkeys(isset_names, False)
        self._is_compiled = False
        
        self.fsf = []
        
        self.log.debug("Starting Group FSF")
        return
    
    def fromDict(self, config):
        self.check_req(config, ["data", "stats"])
        
        # Data
        data = config.pop("data")
        self.check_req(data, ["inputs", "input_type", "model_type", "outdir", "outfsf"])
        ## infile
        self.setData(**data)
        
        # Stats
        stats = config.pop("stats")
        self.check_req(stats, ["evs", "contrasts"])
        self.setStats(**stats)
        
        # Post-Stats
        poststats = config.pop("poststats", {})
        if poststats:
            self.setPostStats(**poststats)
        else:
            self.setPostStats(do_poststats=False)
        
        if len(config) > 0:
            self.log.error("Stuff left-over from config: %s" % config)
        
        return
    
    def run(self, check=True):
        fsf = self.compile()
        self.write(check=check)
        return fsf
    
    def write(self, check=True):
        sfsf = self._remove_extra_stuff("\n".join(self.fsf))
        f = open(self.outfsf, 'w')
        f.write(sfsf)
        f.close()
        if check:
            if 'confoundev_file' in self.fsf_context:
                conf = " " + self.fsf_context['confoundev_file']
            else:
                conf = ""
            self.log.command("feat_model %s%s" % (op.splitext(self.outfsf)[0], conf))
        return
    
    def _remove_extra_stuff(self, a):
        b = re.sub("\n\ +", "\n", a)
        c = re.sub("\n\n+", "\n\n", b)
        return c
    
    def _render(self, fname, **context):
        template = self.env.get_template(fname + '.jinja')
        result = template.render(**context)
        return self._remove_extra_stuff(result)
    
    def compile(self):
        self.log.debug("Compiling fsf...")
        
        # Checks
        self.log.debug("checks")
        for name,isset in self._isset.iteritems():
            if not isset:
                self.log.error("Have not set %s" % name)
        if not len(self.inputs) > 1:
            self.log.fatal("You must specify at least 2 inputs")
        if len(self.evs) == 0:
            self.log.fatal("You must specify at least one EV")
        if len(self.contrasts) == 0:
            self.log.fatal("You must specify at least one contrast")
        
        # Gather fsf        
        fsf = []
        ## header
        self.log.debug("get header")
        fsf_part = self._render('group_header', **self.fsf_context)
        fsf.append(fsf_part)
        ## inputs
        self.log.debug("get inputs")
        fsf_part = self._render('group_inputs', **self.fsf_context)
        fsf.append(fsf_part)
        ## evs
        self.log.debug("get evs")
        fsf_part = self._render('group_evs', **self.fsf_context)
        fsf.append(fsf_part)
        ## cons
        self.log.debug("get cons")
        fsf_part = self._render('group_contrasts', **self.fsf_context)
        fsf.append(fsf_part)
        ## con masks
        self.log.debug("get con masking")
        fsf_part = self._render('group_contrast_masking', n=len(self.contrasts)+len(self.ftests))
        fsf.append(fsf_part)
        ## end
        self.log.debug("get non-gui end")
        fsf_part = self._render('group_nongui', **self.fsf_context)
        fsf.append(fsf_part)
        
        self.fsf = fsf
        return self._remove_extra_stuff("\n".join(fsf))
    
    def setData(self, input_type, inputs, model_type, outdir, outfsf, overwrite=0):
        outdir = self._substitute(outdir)
        inputs = self._getInputsWorker(inputs)
        outfsf = self._substitute(outfsf)
        
        # Input Type
        input_type_choices = {
            'feats': 1,
            'copes': 2
        }
        if input_type not in input_type_choices:
            self.log.error("Input type (%s) must be one of %s" % 
                                (input_type, " ".join(input_type_choices.keys)))
        self.input_type = input_type
        self.fsf_context['input_type'] = input_type_choices[input_type]
        
        # Inputs
        if not len(inputs) > 1:
            self.log.error("Need at least 2 inputs")
        ## check for number of copes and for a reg folder in each input
        if input_type == "feats":
            all_ncopes = [ len(glob(op.join(x, "stats", "cope[0-9]+[.]*"))) for x in inputs ]
            ref_ncopes = all_ncopes[0]
            for i,inpath in enumerate(inputs):
                if not op.isdir(op.join(inpath, "reg")):
                    self.log.warning("Registration directory '%s' does not exist" % 
                                        op.join(inpath, "reg"))
                if all_ncopes[i] != ref_ncopes:
                    self.log.error("Input '%s' has %i copes but should have %i!" % 
                                        (inpath, all_ncopes[i], ref_ncopes))
            self.fsf_context['num_copes'] = ref_ncopes
        self.inputs = inputs
        self.fsf_context['inputs'] = inputs
        self.ninputs = len(inputs)
        self.fsf_context['num_inputs'] = len(inputs)
        
        # Model Type
        model_type_choices = {
            'fixed': 3,
            'ols': 0,
            'flame1': 2,
            'flame2': 1
        }
        if model_type not in model_type_choices:
            self.log.error("Model type (%s) must be one of %s" % 
                                (model_type, " ".join(model_type_choices.keys)))
        self.fsf_context['model_type'] = model_type_choices[model_type]
        
        # Output Directory
        outdir = op.splitext(outdir)[0] + ".gfeat"
        if overwrite == 0 and op.isdir(outdir):
            self.log.error("Output directory '%s' already exists" % outdir)
        elif overwrite > 0:
            self.fsf_context['overwrite'] = overwrite - 1
            overwrite = bool(overwrite - 1)
        if op.isdir(outdir):
            if overwrite:
                self.log.warning("Output directory already exists, will overwrite")
            else:
                self.log.warning("Output directory already exists, instead of overwriting will create a new output directory")
        dir_outdir = op.dirname(outdir)
        if not op.isdir(dir_outdir):
            self.log.info("Creating directory '%s' for feat output" % dir_outdir)
            os.mkdir(dir_outdir)
        self.fsf_context['outputdir'] = outdir
        
        # Output Fsf File
        dir_outfsf = op.dirname(outfsf)
        if not op.isdir(dir_outfsf):
            self.log.info("Creating directory '%s' for fsf output" % dir_outfsf)
            os.mkdir(dir_outfsf)
        self.outfsf = op.splitext(outfsf)[0] + ".fsf"
        
        self._isset['data'] = True
        return
    
    def _getDialect(f):
        sample = f.read(1024)
        f.seek(0)
        
        has_header = csv.Sniffer().has_header(sample)
        if not has_header:
            self.log.warning("file (%s) does not seem to have a header row" % fpath)
        
        dialect = csv.Sniffer().sniff(sample)
        return dialect
    
    def _evsPreCheck(self):
        if not self.inputs:
            self.log.fatal("Must set data before auto-setting EVs")
        if self._isset['evs']:
            self.log.fatal("Can't set EVs more than once")
        return
    
    def _evsPostCheck(self):
        if len(self.evs) == 0:
            self.log.fatal("No EV exist")
        if len(self.groups) == 0:
            self.log.fatal("No group memberships set")
        if len(self.groups) != self.ninputs:
            self.log.error("Should have %i elements for group list but got %i instead" % 
                                (self.ninputs, len(self.groups)))
        for name,vals in self.evs.iteritems():
            if len(vals) != self.ninputs:
                self.log.error("Should have %i elements for EV '%s' but got %i instead." % 
                                    (self.ninputs, name, len(vals)))
        self.fsf_context['groups'] = self.groups
        self.fsf_context['evs'] = self.evs
        self.nevs = len(self.evs)
        self.fsf_context['num_evs'] = self.nevs
        self._isset['evs'] = True
        return
    
    def _setEvsFile(self, evs_file):
        self.log.info("Setting EVs using file '%s'" % evs_file)
        
        # Checks
        evs_file = self._substitute(evs_file)
        if not op.isfile(evs_file):
            self.log.critical("EV file '%s' does not exist" % evs_file)
        self._evsPreCheck()
        
        # Look through file
        with open(evs_file, "rb") as f:
            # Find type of file
            dialect = self._getDialect(f)
            try:
                reader = csv.reader(f, dialect=dialect)
                # Get/check header
                fieldnames = reader.next()
                if fieldnames[0] != "group":
                    self.log.error("First column of EV file '%s' must be 'group'" % evs_file)
                if len(fieldnames) != self.ninputs:
                    self.log.critical("Expected %i columns but got %s in header row" % 
                                        (self.ninputs, len(fieldnames)))
                # Setup
                groups = []
                evs = [ (name, []) for name in fieldnames[1:] ]
                # Go through file
                for row in reader:
                    if len(row) != self.ninputs:
                        self.log.critical("Expected %i columns but got %s (line %d)" % 
                                            (ncols, len(row), reader.line_num))
                    try:
                        row = [ float(r) for r in row ]
                    except ValueError, e:
                        self.log.error("file %s, line %d: something isn't a number" % 
                                            (evs_file, reader.line_num))
                    groups.append(row.pop(0))
                    for ri,r in enumerate(row):
                        evs[ri][1].append(r)
            except csv.Error, e:
                self.log.critical('file %s, line %d: %s' % (evs_file, reader.line_num, e))
            self.groups = groups
            self.evs = OrderedDict(evs)
        self._evsPostCheck()
        return
    
    def _setEvsAuto(self, auto):
        self.log.info("Autosetting EVs (%s)" % auto)
        self._evsPreCheck()
        autoevs = auto.split("_")
        evs = OrderedDict()
        groups = []
        for ev in autoevs:
            if ev == "groupave":
                evs[ev] = [ 1 for x in xrange(self.ninputs) ]
                groups = [ 1 for x in xrange(self.ninputs) ]
            else:
                self.log.error("Unrecognized auto EV '%s'" % ev)
        self.groups = groups
        self.evs = evs
        self._evsPostCheck()
        return
    
    def _contrastsPreCheck(self):
        if not self.evs:
            self.log.fatal("Must set EVs before contrasts")
        if self._isset['contrasts']:
            self.log.fatal("Cannot set contrasts more than once")
    
    def _contrastsPostCheck(self):
        if len(self.contrasts) == 0:
            self.log.fatal("No contrasts exist")
        for name,vals in self.contrasts.iteritems():
            if len(vals) != self.nevs:
                self.log.error("Should have %i elements for contrast '%s' but got %i instead." % 
                                    (self.nevs, name, len(vals)))
        self.ncontrasts = len(self.contrasts)
        self.fsf_context['contrasts'] = self.contrasts
        self.fsf_context['num_cons'] = self.ncontrasts
        self._isset['contrasts'] = True
        return
    
    def _setContrastsFile(self, contrast_file):
        self.log.info("Setting contrasts with file '%s'" % contrast_file)
        
        contrast_file = self._substitute(contrast_file)
        if not op.isfile(contrast_file):
            self.log.critical("Contrast file '%s' does not exist" % contrast_file)
        self._contrastsPreCheck()
        
        with open(evs_file, "rb") as f:
            dialect = self._getDialect(f)
            try:
                reader = csv.reader(f, dialect=dialect)
                contrasts = OrderedDict()
                for row in reader:
                    name = row.pop(0)
                    if len(row) != self.nevs:
                        self.log.error("Expected %i columns for contrast '%s' but got %i" % 
                                            (self.nevs, name, len(row)))
                    try:
                        row = [ float(r) for r in row ]
                    except ValueError, e:
                        self.log.error("file %s, line %d: something isn't a number" % 
                                            (contrast_file, reader.line_num))
                    contrasts[name] = row
            except csv.Error, e:
                self.log.critical('file %s, line %d: %s' % (evs_file, reader.line_num, e))
            self.contrasts = contrasts
        self._contrastsPostCheck()
        return
    
    def _setContrastsAuto(self, auto):
        self.log.info("Autosetting contrasts (%s)" % auto)
        self._contrastsPreCheck()
        autocons = auto.split("_")
        contrasts = OrderedDict()
        for con in autocons:
            if con == "pos":
                contrasts[con] = [1] + [ 0 for x in xrange(self.nevs-1) ]
            elif con == "neg":
                contrasts[con] = [-1] + [ 0 for x in xrange(self.nevs-1) ]
            else:
                self.log.error("Unrecognized auto contrast '%s'" % con)
        self.contrasts = contrasts
        self._contrastsPostCheck()
        return
    
    def setStats(self, evs, contrasts):
        # Evs
        if evs.find("auto_") == 0:
            self._setEvsAuto(evs[5:])
        else:
            self._setEvs(evs)
        # Contrasts
        if contrasts.find("auto_") == 0:
            self._setContrastsAuto(contrasts[5:])
        else:
            self._setContrasts(contrasts)
        # no ftests for now
        self.ftests = []
        self.fsf_context['ftests'] = self.ftests
        self.fsf_context['num_ftests'] = len(self.ftests)
        
        self._isset['stats'] = True
        return
    
    def setPostStats(self, thresh_type="cluster", vox_thresh=2.3, clust_thresh=0.05, prethreshold_mask="", do_poststats=True):
        self.poststats = do_poststats
        
        thresh_type_choices = {
            'none':         0,
            'uncorrected':  1,
            'voxel':        2,
            'cluster':      3
        }
        try:
            self.fsf_context['thresh_type'] = thresh_type_choices[thresh_type]
        except KeyError:
            self.log.error("Unrecognized threshold type '%s'" % thresh_type)
        
        prethreshold_mask = self._substitute(prethreshold_mask)
        self.fsf_context['prethreshold_mask'] = prethreshold_mask
        if prethreshold_mask and not op.isfile(prethreshold_mask):
                self.log.error("Could not find pre-threshold mask %s" % prethreshold_mask)
        
        self._isset['poststats'] = True
        return
    

class FeatGroup(Base):
    """docstring for FeatGroup"""
    
    _logname = "feat_group"
    
    def __init__(self, *args, **kwargs):
        super(FeatGroup, self).__init__(*args, **kwargs)
    
    def fromDict(self, config_dict):
        self.check_req(config_dict, ["infsf"])
        self.setData(**config_dict)
    
    def run(self):
        # Run feat
        self.log.info("Running feat")
        cmd = "feat %s" % self.infsf
        self.log.command(cmd)
        
        # TODO: check input directories
        
        return
    
    def setData(self, infsf):
        infsf = self._substitute(infsf)
        if not op.isfile(infsf):
            self.log.error("Can't find input fsf file '%s'" % infsf)
        self.infsf = infsf
        
        return
    

