import os, re, shutil, tempfile
import os.path as op
from analysis.base import SubjectBase
from collections import OrderedDict
import numpy as np
from subprocess import Popen, PIPE
from copy import deepcopy
from execute import Process

def check_label(label):
    """Display warning if label is greater than 26 characters"""
    l = len(label)
    if l > 24:
        raise DeconWarning(("Label '%s' is too long. " % label) + 
            ("It should have a length of 24 characters or less but has a length of %i." % l))


class CorrelateSubject(SubjectBase):
    
    _logname = "correlate_subject"
    
    def __init__(self, *args, **kwargs):
        super(CorrelateSubject, self).__init__(*args, **kwargs)
        self.log.info("Starting Correlate Subject")
        
        self.infunc = None
        self.outts = None
        self.ts_opts = []
        self._isset_ts = False
        self.ts_cmd = ""
        
        self.cor_opts = []
        self._isset_cor = False
        self.cor_cmd = ""
        return
    
    def fromDict(self, config):
        self.check_req(config, ["cor"])
        
        ts = config.pop("ts", None)
        if ts:
            self.check_req(ts, ["infunc"])
            self.setTs(**ts)
        
        cor = config.pop("cor")
        self.check_req(cor, ["outmap"])
        self.setCor(**cor)
        
        cores = config.pop("cores", 1)
        self.setCores(cores)
        
        for k in config:
            self.log.error("Unknown option '%s' given" % k)
        
        return
    
    def compile(self):
        self.log.info("Compiling")
        if not self._isset_cor:
            self.log.critical("You must set cor before compiling")
        if self.ts_opts:
            self.ts_cmd = " ".join(self.ts_opts)
        self.cor_cmd = " ".join(self.cor_opts)
        return
    
    def run(self):
        self.compile()
        self.log.info("Running")
        
        # 3dmaskave
        if self.ts_cmd:
            f = file(self.outts, 'w')
            self.log.drycommand(self.ts_cmd + " > %s" % self.outts)
            if not self.dry_run:
                p = Process(self.ts_cmd, stdout=f)
                if p.retcode != 0:
                    self.log.error("Error running 3dmaskave")
                    print p.stderr
            f.close()
        
        # 3dTcor1D
        if self.dry_run:
            self.log.drycommand(self.cor_cmd)
        else:
            self.log.command(self.cor_cmd)
        
        return
    
    def setCores(self, val):
        val = str(val)
        self.log.info("Setting cores to %s" % val)
        os.putenv("OMP_NUM_THREADS", val)
        return
    
    def setTs(self, infunc, inmask=None, invox=None, outts=None, **kwargs):
        self.log.info("Setting Time-Series Args/Options")
        ts_opts = ["3dmaskave"]
        
        # Check
        if inmask is not None and invox is not None:
            self.log.error("Cannot specify both input mask (inmask) and input voxel (invox)")
        
        # ROI: Voxel or Mask
        if invox:   # Set voxel coordinates
            if not isinstance(invox, list):
                self.log.error("Input voxel must be a list of 4-5 elements long")
            elif invox[0][1:] == "box":
                if len(invox) != 4:
                    self.log.error("Input box ROI must be a list of 4 elements long")
            elif invox[0][1:] == "ball":
                if len(invox) != 5:
                    self.log.error("Input sphere ROI must be a list of 5 elements long")
            elif invox[0][0] not in ["x", "d", "n", "i"]:
                self.log.error("First element of input voxel list must begin with x, d, n, or i")
            else:
                self.log.error("First element of input voxel list must be Xbox or Xball")
            ts_opts.append("-%s %s" % (invox[0], " ".join(invox[1:])))
        elif inmask: # Or mask
            inmask = self._getInputsWorker(inmask)
            inmask = self._input_list2str(inmask, "input masks")
            if not op.isfile(inmask):
                self.log.error("Input mask '%s' does not exist" % inmask)
            ts_opts.append("-mask %s" % inmask)
        else:       # Eeek
            self.log.error("Must specify either input mask (inmask) or input voxel (invox)")
        
        # Other options
        self.check_ex(kwargs, ["dump", "udump", "indump"])
        for k,v in kwargs.iteritems():
            ts_opts.append("-%s %s" % (k,v))
        
        # Quiet
        ts_opts.append("-quiet")
        
        # Input functional
        infunc = self._getInputsWorker(infunc)
        infunc = self._input_list2str(infunc, "input functionals")
        self.infunc = infunc
        ts_opts.append(infunc)
        
        # Output time-series
        if outts:
            outts = self._substitute(outts)
        else:
            outts = tempfile.mktemp(".1D")  # TODO: when destruct...make sure to delete this
        if not op.isdir(op.dirname(outts)):
            self.log.info("Creating base directory of output '%s'" % op.dirname(outts))
            os.mkdir(op.dirname(outts))
        self.outts = outts
        
        self.ts_opts = ts_opts
        self._isset_ts = True
        return
    
    def setCor(self, outmap, infunc=None, ints=None, mask='', correlation='pearson', outtype='float', overwrite=False):
        self.log.info("Setting Correlation Args/Options")
        cor_opts = ["3dTcorr1D"]
        
        # Options
        ## correlation
        correlation_options = ["pearson", "spearman", "quadrant", "ktaub"]
        if correlation not in correlation_options:
            self.log.error("correlation '%s' not recognized, must be: %s" % 
                                (correlation, ", ".join(correlation_options)))
        cor_opts.append("-%s" % correlation)
        ## output type
        outtype_options = ["float", "short"]
        if outtype not in outtype_options:
            self.log.error("outtype '%s' not recognized, must be: %s" % 
                                (outtype, ", ".join(outtype_options)))
        cor_opts.append("-%s" % outtype)
        ## mask
        mask = self._substitute(mask)
        if mask:
            cor_opts.append("-mask %s" % mask)
        ## output
        outmap = self._substitute(outmap)
        if op.isfile(outmap):
            if overwrite:
                self.log.warning("Removing output '%s'" % outmap)
                os.remove(outmap)
            else:
                self.log.error("Output '%s' already exists" % outmap)
        if not op.isdir(op.dirname(outmap)):
            self.log.info("Creating base directory for output '%s'" % op.dirname(outmap))
            os.mkdir(op.dirname(outmap))
        cor_opts.append("-prefix %s" % outmap)
        
        # Arguments
        ## 4D functional
        ### assign value
        if infunc is None and self.infunc:
            infunc = self.infunc
        elif infunc is None:
            self.log.error("Must specify input functional or must set time-series")
        else:
            infunc = self._getInputsWorker(infunc)
            infunc = self._input_list2str(infunc, "input functionals")
        ### set
        cor_opts.append(infunc)
        ## 1D time-series
        ### assign value
        if ints is None and self.outts:
            ints = self.outts
        elif ints is None:
            self.log.error("Must specify input time-series or must set time-series")
        else:
            ints = self._getInputsWorker(ints)
            ints = self._input_list2str(ints, "input time-series")
            ### check if good
            ncols = self._file_ncols(ints)
            if ncols != 1:
                self.log.error("Input time-series '%s' has %i columns but should have only 1." % 
                                    (ints, ncols))
        ### set
        cor_opts.append(ints)
        
        self.cor_opts = cor_opts
        self._isset_cor = True
        return
    

class DeconSubject(SubjectBase):
    
    _logname = "decon_subject"
    
    def __init__(self, *args, **kwargs):
        super(DeconSubject, self).__init__(*args, **kwargs)
        
        self.evs = OrderedDict()
        self.ev_bases = []
        self.contrasts = OrderedDict()
        self.options = []
        self.cmd_opts = []
        self.no_model = True
        
        isset_names = ['data', 'options', 'ev']
        self._isset = dict.fromkeys(isset_names, False)
    
    def setNoModel(self, val):
        self.no_model = bool(val)
    
    def fromDict(self, config):
        self.check_req(config, ["data", "evs"])
        
        data = config.pop("data")
        self.check_req(data, ["infiles", "outmat", "outpic"])
        self.setData(**data)
        
        motion = config.pop("motion", None)
        if motion:
            self.check_req(motion, ["infiles", "outfile"])
            self.setMotion(**motion)
        
        options = config.pop("options", {})
        self.setOptions(**options)
        
        evs = config.pop("evs")
        self.setEvs(**evs)
        
        contrasts = config.pop("contrasts")
        self.setContrasts(*contrasts)
        
        no_model = config.pop("no_model", None)
        if no_model:
            self.setNoModel(no_model)
        
        for k in config:
            self.log.error("Unknown option '%s' given" % k)
        
        return
    
    def compile(self):
        self.log.info("Compiling")
        
        # Check
        for name,isset in self._isset.iteritems():
            if not isset:
                self.log.error("Have not set %s" % name)
        
        deconvolve_opts = ["3dDeconvolve"]
        
        # Add input
        deconvolve_opts.append("-input %s" % " ".join(self.infiles))
        
        # Add number of EVs
        nevs = len(self.evs)
        deconvolve_opts.append("-num_stimts %s" % nevs)
        
        # Add number of contrasts
        ncons = len(self.contrasts)
        if ncons > 0:
            deconvolve_opts.append("-num_glt %s" % ncons)
        
        # Add options
        deconvolve_opts.extend(self.options)
        
        # Add EVs
        for i,name in enumerate(self.evs):
            ev = self.evs[name]
            ev['i'] = "%02i" % (i+1)
            deconvolve_opts.append("-stim_label %(i)s %(label)s" % ev)
            ev_opt = "-%(option)s %(i)s %(fname)s"
            if ev['model'] is not None:
                ev_opt += " %(model)s"
            deconvolve_opts.append(ev_opt % ev)
        
        # Add contrasts
        if ncons > 0:
            for i,k in enumerate(self.contrasts):
                deconvolve_opts.append("-glt_label %02i %s" % ((i+1), k))
                deconvolve_opts.append("-gltsym 'SYM: %s'" % self.contrasts[k])
        
        # Add stim_bases
        for i in self.ev_bases:
            deconvolve_opts.append("-stim_base %02i" % i)
        
        # Add x1D (matrix)
        deconvolve_opts.append("-x1D %s" % self.outmat)

        # Add xjpeg
        deconvolve_opts.append("-xjpeg %s" % self.outpic)
        
        # x1D_stop
        if self.no_model:
            deconvolve_opts.append("-x1D_stop")
        
        self.cmd_opts = deconvolve_opts
        return " \\\n".join(deconvolve_opts)
    
    def run(self):
        self.compile()
        if self.dry_run:
            self.log.info("Printing command")
            self.log.drycommand(" \\\n".join(self.cmd_opts))
        else:
            self.log.info("Executing command")
            self.log.drycommand(" ".join(self.cmd_opts))
            p = Popen(" ".join(self.cmd_opts), shell=True)  # ghetto fix
            p.communicate()
            if p.returncode != 0:
                self.log.error("Error running 3dDeconvolve")
        return
    
    def setData(self, infiles, outmat, outpic, runfile=None, mkdir=None, overwrite=False):
        self.log.info("Setting Data")
        
        self.infiles = super(DeconSubject, self).getInputs(infiles, runfile, mkdir)
        
        outmat = self._substitute(outmat)
        if op.isfile(outmat):
            if overwrite:
                self.log.warning("Removing previous output matrix: '%s'" % outmat)
                os.remove(outmat)
            else:
                self.log.error("Output matrix '%s' already exists" % outmat)
        self.outmat = outmat
        
        self.outpic = self._substitute(outpic)
        
        self._isset['data'] = True
        return
    
    def setOptions(self, **kwargs):
        self.log.info("Setting Options")
        # todo: check inputs
        for opt,arg in kwargs.iteritems():
            if isinstance(arg, bool) and arg == True:
                self.options.append("-%s" % opt)
            else:
                self.options.append("-%s %s" % (opt, arg))
        
        self._isset['options'] = True
        return
    
    def setEvs(self, stim_file=None, stim_base=None, stim_times=None, **stim_times_kwargs):
        self.log.info("Setting EVs")
                
        if stim_times:
            for x in stim_times:
                name,args = x.items()[0]
                self.addEV("stim_times", name, *args.split())
        
        # Check
        if stim_times_kwargs:
            for k,l in stim_times_kwargs.iteritems():
                if k.find("stim_times_") == -1:
                    self.log.critical("Unrecognized EV option %s" % k)
                for x in l:
                    name,args = x.items()[0]
                    self.addEV(k, name, *args.split())
        
        if stim_file:
            for x in stim_file:
                name,fname = x.items()[0]
                self.addEV("stim_file", name, fname)
        
        if stim_base:
            for name in stim_base:
                self.addStimBase(name)
        
        return
    
    def setContrasts(self, *cons):
        self.log.info("Setting Contrasts")
        for x in cons:
            name,con = x.items()[0]
            self.addContrast(name, con)
        return
    
    def addStimBase(self, name):
        self.log.debug("...adding stim base %s" % name)
        try:
            i = self.evs.keys().index(name) + 1
        except ValueError:
            self.log.error("Could not add -stim_base for %s EV because it doesn't exist" % name)
        
        self.ev_bases.append(i)
        return
    
    def addEV(self, opt, name, fname, model=None):
        self.log.info("...adding EV '%s'" % name)
        # get type
        if opt == "stim_file":
            evtype = 'simple'
        elif opt.find("stim_times") != -1:
            evtype = "complex"
        else:
            self.log.critical("Unrecognized ev option: %s" % opt)
        
        # check name
        check_label(name)
        
        # file
        fname = self._substitute(fname)
        check_fname = re.sub("'?\[.*?\]'?", "", fname)
        if not op.isfile(check_fname):
            self.log.error("Could not find '%s' EV file: %s (%s)" % (name, check_fname, fname))
        
        # model?
        if model is None and evtype == "complex":
            self.log.error("Must specify model for '%s' EV" % name)
        elif model is not None and evtype == "simple":
            self.log.error("Cannot specify model for '%s' EV" % name)
        elif model is not None and model.find("'") != 0:
            model = "'%s'" % model
        
        # save
        if name in self.evs:
            self.log.error("'%s' EV has been specified more than once" % name)
        self.evs[name] = {
            'option': opt,
            'label': name,
            'fname': fname,
            'model': model
        }
        
        self._isset['ev'] = True
        return
    
    def addContrast(self, name, con):
        self.log.info("...adding contrast %s" % name)
        check_label(name)
        if name in self.contrasts:
            self.log.error("%s contrast has been specified more than once" % name)
        self.contrasts[name] = con
        
        self._isset['contrast'] = True
        return
    


class RemlSubject(SubjectBase):
    
    _logname = "reml_subject"
    
    def __init__(self, *args, **kwargs):
        super(RemlSubject, self).__init__(*args, **kwargs)
        
        self.options = []
        self.cmd_opts = []
        self.type = 'reml'
        
        isset_names = ['data']
        self._isset = dict.fromkeys(isset_names, False)
    
    def setType(self, val):
        self.log.info("Setting Output Type %s" % val)
        if val != "reml" and val != "ols":
            self.log.error("Incorrect REML type '%s'. Must be 'reml' or 'ols'" % val)
        self.type = val
    
    def setCores(self, val):
        val = str(val)
        self.log.info("Setting cores to %s" % val)
        os.putenv("OMP_NUM_THREADS", val)
        return
    
    def fromDict(self, config):
        self.check_req(config, ["data"])
        
        data = config.pop("data")
        self.check_req(data, ["infiles", "inmat", "outdir"])
        self.setData(**data)
        
        options = config.pop("options", {})
        self.setOptions(**options)
        
        otype = config.pop("type", None)
        if otype:
            self.setType(otype)
        
        cores = config.pop("cores", 1)
        self.setCores(cores)
        
        for k in config:
            self.log.error("Unknown option '%s' given" % k)
        
        return
    
    def compile(self):
        self.log.info("Compiling")
        
        # Check
        for name,isset in self._isset.iteritems():
            if not isset:
                self.log.error("Have not set %s" % name)
        
        reml_opts = ["3dREMLfit",]

        reml_opts.append("-input '%s'" % " ".join(self.infiles))
        reml_opts.append("-matrix %s" % self.inmat)

        if self.mask:
           reml_opts.append("-mask %s" % self.mask)

        # What type of outputs
        if self.type == "reml":
           outputs = [
               "-Rbuck stats/bucket.nii.gz",
               "-Rbeta stats/betas.nii.gz",
               "-Rvar stats/variances.nii.gz",
               "-Rfitts stats/fitts.nii.gz",
               "-Rerrts stats/residuals.nii.gz",
               "-Rwherr stats/whitened_residuals.nii.gz",
           ]
        elif self.type == "ols":
           outputs = [
               "-Obuck stats/bucket.nii.gz",
               "-Obeta stats/betas.nii.gz",
               "-Ofitts stats/fitts.nii.gz",
               "-Oerrts stats/residuals.nii.gz"
           ]

        reml_opts.extend(outputs)
        self.cmd_opts = reml_opts

        return "\\\n".join(reml_opts)

    def setupOutput(self):
        """Sets up the contents of subject output directory. ds => dictionary"""
        self.log.info("Setting up contents of output directory")
        
        dirs = ["custom_timing_files", "stats", "reg_standard", "tsplot", "statplot", "logs"]
        
        self.log.info("...creating output directory")
        os.mkdir(self.outdir)
        
        self.log.info("...copying input matrix")
        shutil.copy(self.inmat, op.join(self.outdir, "model.1D"))
        
        # inputs linked
        self.log.info("...soft linking inputs")
        os.mkdir(op.join(self.outdir, "inputs"))
        for i in range(len(self.infiles)):
            ext = op.splitext(self.infiles[i])[1]
            if ext == ".gz":
                ext = op.splitext(op.splitext(self.infiles[i])[0])[1] + ext
            os.symlink(self.infiles[i], 
                op.join(self.outdir, "inputs", "filtered_func_run%02i%s" % (i+1, ext)))
                
        # copy mask
        if self.mask:
            if isinstance(self.mask, str):
                self.log.info("...copying mask")
                ext = op.splitext(self.mask)[1]
                if ext == ".gz":
                    ext = op.splitext(op.splitext(self.mask)[0])[1] + ext
                shutil.copy(self.mask, op.join(self.outdir, "mask%s" % ext))
            elif isinstance(self.mask, list):
                ext = op.splitext(self.mask[0])[1]
                if ext == ".gz":
                    ext = op.splitext(op.splitext(self.mask[0])[0])[1] + ext
                self.log.info("...merging masks across subjects")
                self.log.command("fslmerge -t %s %s" % (
                    op.join(self.outdir, "allsubjs_masks%s" % ext),
                    " ".join(self.mask)
                ))
                self.log.command("fslmaths %s -Tmin -bin %s" % (
                    op.join(self.outdir, "allsubjs_masks%s" % ext),
                    op.join(self.outdir, "mask%s" % ext),
                ))
            self.mask = op.join(self.outdir, "mask%s" % ext)
        
        # -> custom_timing_files
        # -> stats folder
        # -> reg_standard
        # -> tsplot
        # -> statplot
        self.log.info("...creating various directories")
        for d in dirs:
            os.mkdir(op.join(self.outdir, d))
        
        return
    
    def run(self):
        # setup output
        self.setupOutput()
        # compile
        self.compile()
        # run
        if self.dry_run:
            self.log.drycommand("\\\n".join(self.cmd_opts))
        else:
            self.log.info("Running stats")
            self.log.drycommand(" ".join(self.cmd_opts))
            p = Popen(" ".join(self.cmd_opts), shell=True, cwd=self.outdir)  # ghetto fix
            p.communicate()
            if p.returncode != 0:
                self.log.error("Error running 3dREMLfit")
            f = file(op.join(self.outdir, "logs", "reml.cmd"), 'w')
            f.write("\\\n".join(self.cmd_opts))
            f.close()
    
    def setData(self, infiles, inmat, outdir, mask=None, runfile=None, mkdir=None, overwrite=0):
        self.log.info("Setting Data")
        
        self.infiles = super(RemlSubject, self).getInputs(infiles, runfile, mkdir)
        
        self.inmat = self._substitute(inmat)
        if not op.isfile(self.inmat):
            self.log.error("Input matrix from 3dDeconvolve '%s' not found." % self.inmat)
        
        outdir = self._substitute(outdir)
        outdir = "%s.reml" % (op.splitext(outdir)[0])
        if op.isdir(outdir):
            if overwrite == 0:
                self.log.error("Output '%s' already exists" % outdir)
            elif overwrite == 1:
                outdir = "%s+.reml" % (op.splitext(outdir)[0])
                self.log.warning("Output already exists. Saving as '%s'." % outdir)
            else:
                self.log.warning("Removing output directory '%s'." % outdir)
                shutil.rmtree(outdir)
        self.outdir = outdir
        
        if mask:
            mask = self._getInputsWorker(mask, 'file')
        self.mask = mask
        
        self._isset['data'] = True
        return
    
    def setOptions(self, **kwargs):
        self.log.info("Setting Options")
        # todo: check inputs
        for opt,arg in kwargs.iteritems():
            if isinstance(arg, bool) and arg == True:
                self.options.append("-%s" % opt)
            else:
                self.options.append("-%s %s" % (opt, arg))
        
        self._isset['options'] = True
        return
    

# class for beta-series correlations?    
class BetaSeriesSubject(SubjectBase):
    
    _logname = "betaseries_subject"
    
    def __init__(self, *args, **kwargs):
        super(BetaSeriesSubject, self).__init__(*args, **kwargs)
        self.log.info("Starting Creation of Beta-Series")
        
        self.evs = OrderedDict()
        isset_names = ['data', 'evs']
        self._isset = dict.fromkeys(isset_names, False)
    
    def fromDict(self, config):
        self.check_req(config, ["data", "evs"])
        
        data = config.pop("data")
        self.check_req(data, ["indir"])
        self.setData(**data)
        
        evs = config.pop("evs")
        self.setEvs(*evs)
        
        for k in config:
            self.log.error("Unknown option '%s' given" % k)
        
        return
    
    def compile(self):
        self.log.info("Compiling")
        
        # Check
        for name,isset in self._isset.iteritems():
            if not isset:
                self.log.error("Have not set %s" % name)
        
        return
    
    def run(self):
        self.compile()
        # Save beta-series for EVs
        p = self.log.command("3dAttribute BRICK_LABS stats/betas.nii.gz", cwd=self.indir)
        all_labels = p.stdout.split("~")[0:-1]
        evi = 0
        for name,step in self.evs.iteritems():
            evi += 1
            label_inds = np.array([ i for i,x in enumerate(all_labels) if x.find(name) == 0 ])
            for start in xrange(step):
                filt_label_inds = [ str(x) for x in label_inds[range(start, label_inds.shape[0], step)] ]
                fname = "ev%02i_%s_num%i.nii.gz" % (evi, name, start+1)
                self.log.command("3dcalc -a stats/betas.nii.gz'[%s]' -expr a -prefix betaseries/%s" % (",".join(filt_label_inds), fname), cwd=self.indir)
                self.log.command("buc2func.R betaseries/%s betaseries/%s" % (fname, fname), cwd=self.indir)
        
        # applywarp -i ${rundir}/${inf} -r ${regdir}/standard.nii.gz -o ${rundir}/${outf} -w ${regdir}/highres2standard_warp.nii.gz --premat=${regdir}/example_func2highres.mat --interp=spline -v
    def setData(self, indir, regdir=""):
        self.log.info("Setting data")
        
        indir = self._substitute(indir)
        regdir = self._substitute(regdir)
        
        if not op.isdir(indir):
            self.log.error("Directory '%s' does not exist" % indir)
        if not op.isdir(op.join(indir, "betaseries")):
            self.log.info("...creating betaseries directory")
            os.mkdir(op.join(indir, "betaseries"))
        if regdir and not op.isdir(regdir):
            self.log.error("Reg directory '%s' does not exist" % regdir)
            
        self.indir = indir
        self.regdir = regdir
        
        self._isset['data'] = True
        return
    
    def addEV(self, name, steps):
        if name in self.evs:
            self.log.error("Duplicate EV %s" % name)
        self.evs[name] = steps
        
        self._isset['evs'] = True
        return
    
    def setEvs(self, *evs):
        self.log.info("Setting EVs")
        for ev in evs:
            name,step = ev.items()[0]
            self.addEV(name, step)
        return
    

class RegisterBetaSeriesSubject(SubjectBase):
    
    _logname = "registerbetaseries_subject"
    _std = "standard.nii.gz"
    _warp = "highres2standard_warp.nii.gz"
    _mat = "example_func2highres.mat"
    
    def __init__(self, *args, **kwargs):
        super(RegisterBetaSeriesSubject, self).__init__(*args, **kwargs)
        self.log.info("Starting Creation of Beta-Series")        
        isset_names = ['data', 'options']
        self._isset = dict.fromkeys(isset_names, False)
    
    def fromDict(self, config):
        self.log.debug("beginning setup")
        
        self.check_req(config, ["data"])
        
        data = config.pop("data")
        self.check_req(data, ["indir", "regdir"])
        self.setData(**data)
        
        options = config.pop("options", {})
        self.setOptions(**options)
        
        for k in config:
            self.log.error("Unknown option '%s' given" % k)
        
        return
    
    def compile(self):
        self.log.info("Compiling")
        # Check
        for name,isset in self._isset.iteritems():
            if not isset:
                self.log.error("Have not set %s" % name)
        # Commands
        commands = []
        if self.bsdir == self.outdir:
            raise Exception("todo")
        for infname in self.bsfnames:
            cmd = ["applywarp"]
            cmd.append("-i %s" % op.join(self.bsdir, infname))
            cmd.append("--premat=%s" % op.join(self.regdir, self._mat))
            cmd.append("-w %s" % op.join(self.regdir, self._warp))
            cmd.append("-r %s" % op.join(self.regdir, self._std))
            cmd.append("-o %s" % op.join(self.outdir, infname))
            cmd.extend(self.cmd_opts)
            commands.append(" ".join(cmd)) 
        self.commands = commands
        return       
    
    def run(self):
        self.compile()
        self.log.info("Running")
        for cmd in self.commands:
            if self.dry_run:
                self.log.drycommand(cmd)
            else:
                self.log.command(cmd)
        return
    
    def setData(self, indir, regdir, overwrite=False):
        self.log.info("Setting data")
        
        indir = self._substitute(indir)
        regdir = self._substitute(regdir)
        bsdir = op.join(indir, "betaseries")
        outdir = op.join(indir, "betaseries2standard")
        
        
        if not op.isdir(indir):
            self.log.error("Directory '%s' does not exist" % indir)
        if not op.isdir(regdir):
            self.log.error("Reg directory '%s' does not exist" % regdir)
        if not op.isdir(bsdir) or len(os.listdir(bsdir)) == 0:
            self.log.error("You must first extract out the different beta-series")
        if op.isdir(outdir) and len(os.listdir(outdir)) > 0:
            if overwrite:
                self.log.warning("Removing registration of beta-series directory")
                shutil.rmtree(outdir)
            else:
                self.log.error("Registration of beta-series has been done")
        elif not op.isdir(outdir):
            os.mkdir(outdir)
        
        self.indir = indir
        self.regdir = regdir
        self.bsdir = bsdir
        self.outdir = outdir
        
        bsfnames = os.listdir(self.bsdir)
        bsfnames.sort()
        self.bsfnames = bsfnames
        
        self._isset['data'] = True
        return
    
    def setOptions(self, **kwargs):
        self.log.info("Setting Options")
        
        cmd_opts = []
        for k,v in kwargs.iteritems():
            if len(k) == 1:
                pre = "-"
            else:
                pre = "--"
            if isinstance(v, bool) and v == True:
                cmd_opts.append("%s%s" % (pre, k))
            else:
                cmd_opts.append("%s%s=%s" % (pre, k, v))
        self.cmd_opts = cmd_opts
        
        self._isset['options'] = True
        return
    
    

# # -> reg folder linked
# self.log.info("...soft linking reg directory")
# os.symlink(self.regdir, op.join(self.outdir, "reg"))