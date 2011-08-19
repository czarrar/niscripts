# TODO: split this stuff up

import os, nibabel, re, yaml, tempfile
import os.path as op
from string import Template
from collections import OrderedDict
from jinja2 import Environment, PackageLoader
from glob import glob
import numpy as np
from copy import deepcopy
from execute import Process
from analysis.base import SubjectBase

# Time
NIFTI_UNITS_SEC=8
NIFTI_UNITS_MSEC=16
NIFTI_UNITS_USEC=24

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
        
        decon = config.pop("decon", {})
        self.setDecon(**decon)
        
        return
    
    def compile(self):
        pass
    
    def run(self):
        if not self._isset_data or not self._isset_decon:
            self.log.fatal("Data or 3dDeconvolve must be setup before running")
        
        # Combine Functionals
        cmd = "fslmerge -t %s %s" % (self.outfunc, " ".join(self.infiles))
        self.log.command(cmd, cwd=op.dirname(self.outfunc))
        
        # 3dDeconvolve
        if self.decon_opts:
            if not self.outconfound:
                self.log.error("Must specify outconfound file if doing 3dDeconvole")
            # deconvolve
            self.decon_opts.insert(0, "-input %s" % " ".join(self.infiles))
            self.decon_opts.insert(0, "3dDeconvolve")
            self.log.command(" ".join(self.decon_opts), cwd=op.dirname(self.decon_outmat))
            # confound file
            x = np.loadtxt(self.decon_outmat)
            if self.inconfound:
                x2 = np.loadtxt(self.inconfound)
                if x.shape[0] != x2.shape[0]:
                    self.log.error("inconfound '%s' and 3dDeconvolve output '%s' do not have same number of rows" % (self.inconfound, self.decon_outmat))
                x3 = np.zeros((x.shape[0], x.shape[1]+x2.shape[1]))
                x3[:,0:x.shape[1]] = x
                x3[:,x.shape[1]:] = x2
                x = x3
            np.savetxt(self.outconfound, x, fmt="%.9f")
            os.remove(self.decon_outmat)
        
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
    
    def setDecon(self, polort=None, inconfound="", outconfound="", tr=None):
        """docstring for setDecon"""
        self.log.debug("Setting 3dDeconvolve")
        decon_opts = []     # add -input later
        
        inconfound = self._substitute(inconfound)
        outconfound = self._substitute(outconfound)
        
        if polort:
            decon_opts.append("-polort %s" % polort)
        if decon_opts:
            if tr:
                decon_opts.append("-force_TR %f" % float(tr))
            outmat = tempfile.mktemp(".1D")
            decon_opts.append("-x1D %s" % outmat)
            self.decon_outmat = outmat
            if inconfound and not op.isfile(inconfound):
                self.log.error("Cannot find inconfound file '%s'" % inconfound)
            if not outconfound:
                self.log.error("You need to specify outconfound file to run 3dDeconvolve")
            decon_opts.append("-x1D_stop")
        elif inconfound:
            self.log.error("Cannot specify inconfound without polort or censortr option")
        
        self.inconfound = inconfound
        self.outconfound = outconfound
        self.decon_opts = decon_opts
        self.tr = tr
        
        self._isset_decon = True
        
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
        ## infile
        infile = data.pop("infile")
        infile = self._substitute(infile)
        gpaths = glob(infile)
        if len(gpaths) == 1:
            infile = gpaths[0]
        elif len(gpaths) > 1:
            self.log.error("Cannot have more than one infile %s from %s" % (gpaths, infile))
        else:
            self.log.error("Did not find infile %s" % infile)
        ## outdir
        outdir = data.pop("outdir")
        outdir = self._substitute(outdir)
        self.setData(infile, outdir, **data)
        
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
    
    def setData(self, infile, outdir, outfsf, highpass_filter, tr=None, overwrite=0):
        infile = self._substitute(infile)
        outdir = self._substitute(outdir)
        outfsf = self._substitute(outfsf)
        
        if not op.isfile(infile):
            self.log.error("Could not find input file %s" % infile)
        
        self.fsf_context['func_file'] = infile
        self.fsf_context['outputdir'] = outdir
        self.fsf_context['highpass_yn'] = False
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
    
    def setStats(self, prewhiten=True, confoundev_file="", evs={}, contrasts={}):
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
    
    
    def addEV(self, name, fname, convolution, tempfilt=False, tempderiv=False, **opts):
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
        
        # Waveform Shape (1 column or 3 column file)
        ncols = self._file_ncols(fname)
        if ncols == 1:
            waveform = waveform_choices['custom1']
        elif ncols == 3:
            waveform = waveform_choices['custom3']
        if ncols != 1 and ncols != 3:
            self.log.error("File %s for EV %s must have 1 or 3 columns, instead %i were found" % 
                                (fname, name, ncols))
        
        # Convolution
        try:
            convolution = convolution_choices[convolution]
        except KeyError as k:
            self.log.error("Convolution '%s' not recognized" % convolution)
        
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
    

class FsfInfo(object):
    def __init__(self, fname):
        """docstring for __init__"""
        f = file(fname, 'r')
        self.fsf = [ l.strip() for l in f if l.strip() != '' and l[0] != '#' ]
        return
    
    def getOutput(self):
        res = [ re.search("set\ fmri\(.*?\)\ (.*)", l) for l in fsf if l.find("outputdir") != -1 ]
        return res.groups()[0].strip('"')

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
        outputdir = fsfinfo.getOutput()
        
        # Soft link regdir
        self.log.info("Soft linking registration directory")
        os.symlink(self.regdir, op.join("outputdir", "reg"))
        
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

class RegressSubject(SubjectBase):
    
    _logname = "regress_subject"
    
    def __init__(self, *args, **kwargs):
        super(FeatSubject, self).__init__(*args, **kwargs)
    
    def fromDict(self, config_dict):
        self.check_req(config_dict, ["design_file", "in_file", "out_file"])
        self.cmd_opts = self.setVars(config_dict)
    
    def setVars(self, config_dict):
        choices = ["in_file", "out_file", "mask"]
        for k in choices:
            if k in config_dict:
                config_dict[k] = self._substitute(config_dict[k])
        return config_dict
    
    def run(self):
        self.log.info("Setting up command")
        filt = fsl.FilterRegressor(**cmd_opts)
        cmd = filt.cmdline
        
        self.log.info("Running fsl_regfilt")
        self.log.command(cmd)
    

