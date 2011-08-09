#!/usr/bin/env python

import argparse, os, sys, yaml, shutil
import os.path as op
from copy import deepcopy
from glob import glob
from string import Template
import zlogger
from usage import NiArgumentParser

#####
# Process user arguments
#####

def create_parser():
    parser = NiArgumentParser(fromfile_prefix_chars='@', 
                description="GLM/REML analysis of subject level activations")
    
    group = parser.add_argument_group('Required')
    group.add_argument('-s', '--subjects', nargs="+", required=True)
    group.add_argument('-c', '--config', type=argparse.FileType('r'), required=True)
    
    group = parser.add_argument_group('Optional')
    group.add_argument("--interp", choices=["lin", "nn", "sinc", "spline"], default="lin")
    group.add_argument("--dry-run", action="store_true", default=False)
    group.add_argument("--cores", type=int, default=1)
    group.add_argument("--verbose", action="store_const", const=1, dest="verbosity", default=0)
    group.add_argument("--debug", action="store_const", const=2, dest="verbosity", default=0)
    group.add_argument("--force", action="store_true", default=False, help="doesn't work")
    
    return parser


def test_wf(dry_run=True):
    if dry_run:
        arglist = "-c %s/tests/deconv.yaml --dry-run --cores 1 --force --debug -s tb3417" % os.getenv("NISCRIPTS")
    else:
        arglist = "-c %s/tests/deconv.yaml --cores 1 --debug --force -s tb3417" % os.getenv("NISCRIPTS")
    return main(arglist.split())

def main(arglist):
    # Parse
    parser = create_parser()
    args = parser.parse_args(arglist)
    workflow = Analyze(args.subjects, args.config, args.interp, args.cores, args.dry_run, args.verbosity, args.force)
    workflow.run()
    return workflow

#####
# Needed Classes
#####

class DeconError(Exception):
    """General error when running this script"""

class DeconWarning(Exception):
    """General warning when running this script"""

class Struct:
    def __init__(self, **entries):
        """Can easily convert a dictionary to an object
        
        >> d = {'one': 1, 'two': 2}
        >> obj = Struct(**d)
        >> print(obj.one)
        >> print(obj.two)
        
        """
        self.__dict__.update(entries)
    

class AnalyzeSubject(object):
    model_matrix = "model.1D"
    model_matrix_pic = "model.jpg"
    
    def __init__(self, vars, decon, reml, loglevel=zlogger.logging.DEBUG, interp="lin", dry_run=False):
        # todo: check vars has subject, func...
        self.__dict__.update(**vars)
        self.vars = vars
        # create output
        if dry_run is True:
            log_fname = None
        else:
            log_fname = op.join(self.outputdir, "logs", "glm_subject.log")
            os.mkdir(self.outputdir)
            os.mkdir(op.join(self.outputdir, "logs"))
        # create logger
        self.slog = zlogger.getLogger("glm_%s" % self.subject, level=loglevel, 
                                      fname=log_fname)
        # save args
        self.tdecon = deepcopy(decon)    # template
        self.treml = deepcopy(reml)      # template
        self.interp = interp
        self.dry_run = dry_run
        return
    
    def __call__(self):
        return self.run()
    
    def run(self):
        print
        
        if not self.dry_run:
            self.setup_output()
            print
        
        # 3dDeconvolve
        decon = self.create_model()
        if self.dry_run is True:
            self.slog.drycommand(" \\\n".join(decon))
        else:
            p = self.slog.command(" ".join(decon), cwd=self.outputdir)
            if p.retcode != 0:
                raise DeconError("Non-zero exit from 3dDeconvolve")
        print
        
        # 3dREMLfit
        reml = self.run_model()
        if self.dry_run is True:
            self.slog.drycommand(" \\\n".join(reml))
        else:
            p = self.slog.command(" ".join(reml), cwd=self.outputdir)
            if p.retcode != 0:
                raise DeconError("Non-zero exit from 3dREMLfit")
        print
        
        # Extract GLTs
        self.extract_glts()
        print
        
        # Apply Registration
        self.apply_reg()
        print
        
        return
    
    def setup_output(self):
        """Sets up the contents of subject output directory. ds => dictionary"""
        self.slog.title("Setting up contents of output directory")
        dirs = ["custom_timing_files", "stats", "reg_standard", "tsplot", "statplot"]
        
        # inputs linked
        self.slog.info("...soft linking inputs")
        os.mkdir(op.join(self.outputdir, "inputs"))
        for i in range(len(self.func)):
            ext = op.splitext(self.func[i])[1]
            os.symlink(self.func[i], 
                op.join(self.outputdir, "inputs", "filtered_func_run%02i%s" % (i+1, ext)))
        
        # -> reg folder linked
        self.slog.info("...soft linking reg directory")
        os.symlink(self.regdir, op.join(self.outputdir, "reg"))
        
        # copy mask
        ext = op.splitext(self.mask)[1]
        if ext == "gz":
            ext = op.splitext(op.splitext(self.mask)[0])[1] + ext
        if isinstance(self.mask, str):
            self.slog.info("...copying mask")
            shutil.copy(self.mask, op.join(self.outputdir, "mask%s" % ext))
        elif isinstance(self.mask, list):
            self.slog.info("...merging masks across subjects")
            self.slog.command("fslmerge -t %s %s" % (
                op.join(self.outputdir, "allsubjs_masks%s" % ext),
                " ".join(self.mask)
            ))
            self.slog.command("fslmaths %s -Tmin -bin %s" % (
                op.join(self.outputdir, "allsubjs_masks%s" % ext),
                op.join(self.outputdir, "mask%s" % ext),
            ))
        
        # -> custom_timing_files
        # -> stats folder
        # -> reg_standard
        # -> tsplot
        # -> statplot
        self.slog.info("...creating various directories")
        for d in dirs:
            os.mkdir(op.join(self.outputdir, d))
        
        return
    
    def create_model(self):
        """Generate command to run 3dDeconvolve"""
        
        self.slog.info("Creating model")
        
        # gather options
        deconvolve_opts = []
        deconvolve = self.tdecon
        check_ex(deconvolve, [
            "input", "num_stimts", "num_glt", "noFDR", 
            "x1D", "xjpeg", "x1D_stop", "nox1D"
        ])
        
        # check for polort, force_TR, local_times
        copts = ["polort", "force_TR", "local_times"]
        for k in copts:
            if k in deconvolve:
                v = deconvolve[k]
                if isinstance(v, str) and v == ".":
                    deconvolve_opts.append("-%s" % k)
                else:
                    deconvolve_opts.append("-%s %s" % (k, v))
        
        # Go through stim_times, stim_times_ETC, stim_file, 
        # TODO: record any stim files and copy them over
        stim_labels = []
        ks = [ x for x in deconvolve.keys() if x.find("stim_times") == 0 ]
        if "stim_file" in deconvolve:
            ks.append("stim_file")
        for k in ks:
            opts = deconvolve.pop(k)
            for e in opts:
                if len(e) > 1:
                    pass # raise error
                label, args = e.items()[0]
                check_label(label)
                try:
                    args = zsub_basic(args, self.vars)
                except KeyError as strerr:
                    raise DeconError("Could not find the variable %s in %s with label %s" % (
                        strerr, args, k))
                stim_labels.append(label)
                si = len(stim_labels)
                # e.g. -stim_label 01 name
                deconvolve_opts.append("-stim_label %02i %s" % (si, label))
                # e.g. -stim_times 01 fname model
                deconvolve_opts.append("-%s %02i %s" % (k, si, args))

        # Go through stim_base
        if "stim_base" in deconvolve:
            labels = deconvolve.pop("stim_base")
            if isinstance(labels, str):
                opts = [labels]
            elif not isinstance(opts, list):
                pass # raise an error
            for label in labels:
                check_label(label)
                try:
                    si = stim_labels.index(label) + 1
                    deconvolve_opts.append("-stim_base %02i" % si)
                except ValueError:
                    raise DeconError("todo")
                
        # Go through gltsym
        glt_labels = []
        if "gltsym" in deconvolve:
            opts = deconvolve.pop("gltsym")
            for e in opts:
                if len(e) > 1:
                    pass # raise error
                if not isinstance(e, dict):
                    raise DeconError("Line with gltsym in config file is not properly formatted: %s" % e)
                label, args = e.items()[0]
                check_label(label)
                glt_labels.append(label)
                si = len(glt_labels)
                # e.g., -glt_label 01 Memoranda
                deconvolve_opts.append("-glt_label %02i %s" % (si, label))
                # e.g., -gltsym 'SYM: +1*Memoranda'
                deconvolve_opts.append("-gltsym 'SYM: %s'" % args)
        
        # Add number of glts
        self.glt_labels = glt_labels
        if len(glt_labels) > 0:
            deconvolve_opts.insert(0, "-num_glt %s" % len(glt_labels))
        
        # Add number of stims
        self.stim_labels = stim_labels
        if len(stim_labels) > 0:
            deconvolve_opts.insert(0, "-num_stimts %s" % len(stim_labels))
        else:
            raise DeconError("You must specify some evs/stimuli")
        
        # Add input
        deconvolve_opts.insert(0, "-input %s" % " ".join(self.func))
        
        # Add noFDR
        deconvolve_opts.append("-noFDR")
        
        # Go through everything else
        ## TODO: add checks that any other option is truly a 3dDeconvolve option
        for k,v in deconvolve.iteritems():
            if isinstance(v, str) and v == ".":
                deconvolve_opts.append("-%s" % k)
            else:
                deconvolve_opts.append("-%s %s" % (k, v))
        
        # Add x1D (matrix)
        deconvolve_opts.append("-x1D %s" % op.join(self.outputdir, self.model_matrix))

        # Add xjpeg
        deconvolve_opts.append("-xjpeg %s" % op.join(self.outputdir, self.model_matrix_pic))
        
        # x1D_stop
        deconvolve_opts.append("-x1D_stop")
        
        deconvolve_opts.insert(0, "3dDeconvolve")
        
        return deconvolve_opts
    
    def run_model(self):
        """Generate command to run 3dREMLfit (depends on stuff from deconvolve)"""
        
        self.slog.info("Running model")
        
        reml_opts = []
        reml = self.treml
        # todo: check_ex(reml, ["Rvar"...etc])
        
        reml_opts.append("-input '%s'" % " ".join(self.func))
        reml_opts.append("-matrix %s" % op.join(self.outputdir, self.model_matrix))
        reml_opts.append("-mask %s" % self.mask)
        
        # What type of outputs
        if "type" in reml:
            glm_type = reml.pop("type")
            if glm_type == "reml":
                outputs = [
                    "-Rbuck stats/bucket.nii.gz",
                    "-Rbeta stats/betas.nii.gz",
                    "-Rvar stats/variances.nii.gz",
                    "-Rfitts stats/fitts.nii.gz",
                    "-Rerrts stats/residuals.nii.gz",
                    "-Rwherr stats/whitened_residuals.nii.gz",
                ]
            elif glm_type == "ols":
                outputs = [
                    "-Obuck stats/bucket.nii.gz",
                    "-Obeta stats/betas.nii.gz",
                    "-Ofitts stats/fitts.nii.gz",
                    "-Oerrts stats/residuals.nii.gz"
                ]
            else:
                raise DeconError("run_model: type can only have 'reml' or 'ols' as values")
            reml_opts.extend(outputs)
            self.glm_type = glm_type
        else:
            self.glm_type = None
        
        reml_opts.extend([
            "-tout",
            "-noFDR",
            "-verb"
        ])
        
        reml_opts.insert(0, "3dREMLfit")
        
        return reml_opts
    
    def extract_glts(self):
        # 1. Get glts bucket
        if not self.dry_run:
            self.slog.info("Extracting all GLTs into seperate bucket")
            p = self.slog.command("3dinfo -label2index '%s#0_Coef' %s" % 
                                            (self.glt_labels[0], "stats/bucket.nii.gz"), 
                                        cwd=self.outputdir)
            if p.retcode != 0:
                raise DeconError("Error extracting label index")
            i = int(p.stdout)
            self.slog.command("3dcalc -a %s'[%i..$]' -expr 'a' -prefix %s" % 
                                ("stats/bucket.nii.gz", i, "stats/glts.nii.gz"), 
                             cwd=self.outputdir)
        
        # 2. Get individual glt coefs and tstats
        self.slog.info("Extracting individual GLT coefs and tstats")
        for i in xrange(len(self.glt_labels)):
            label = self.glt_labels[i]
            ## coef
            cmd = "3dcalc -a %s'[%s#0_Coef]' -expr 'a' -prefix stats/coef%02i_%s.nii.gz" % ("stats/bucket.nii.gz", label, i, label)
            if self.dry_run:
                self.slog.drycommand(cmd)
            else:
                p = self.slog.command(cmd, cwd=self.outputdir)
                if p.retcode != 0:
                    raise DeconError("Error extracting Coef for %s from bucket" % label)
            ## tstat
            cmd = "3dcalc -a %s'[%s#0_Tstat]' -expr 'a' -prefix stats/tstat%02i_%s.nii.gz" % ("stats/bucket.nii.gz", label, i, label)
            if self.dry_run:
                self.slog.drycommand(cmd)
            else:
                p = self.slog.command(cmd, cwd=self.outputdir)
                if p.retcode != 0:
                    raise DeconError("Error extracting Tstat for %s from bucket" % label)
    
    def easy_thresh(self):
        """For now this is not applied"""
        # 1. smoothest
        ## dof
        p = Process("3dAttribute BRICK_STATAUX %s'[0]'" % "stats/bucket.nii.gz", cwd=self.outputdir) | Process("awk '{print $5}'", cwd=self.outputdir)
        if p.retcode != 0:
            raise DeconError("couldn't get DOF from bucket")
        dof = p.stdout
        ## res4d
        if self.glm_type == "reml":
            res4d = "stats/whitened_residuals.nii.gz"
        elif self.glm_type == "ols":
            rest4d = "stats/residuals.nii.gz"
        else:
            raise Exception("todo")
        ## run
        p = self.slog.command("smoothest -d %s -m %s -r %s > stats/smoothness" % (dof, self.mask, res4d), cwd=self.outputdir)
        if p.retcode != 0:
            raise DeconError("failure running FSL's smoothest")
        dlh, volume, resels = tuple([ x.split(" ")[1] for x in p.stdout.splitlines() ])
        
        # 2. apply cluster
        #cluster -i ${threshZstat} -c ${cope} -t ${voxThresh} -p ${clustThresh} -d \${DLH} --volume=\${VOLUME} --othresh=${threshZstat} -o ${threshDir}/cluster_mask_zstat${i} --olmax=${threshDir}/lmax_zstat${i}.txt > ${threshDir}/cluster_zstat${i}.txt
    
    def apply_reg(self):
        self.slog.info("Applying registration to Coefs and Tstats")
        
        # fnirt or flirt?
        if op.isfile(op.join(self.outputdir, "reg", "highres2standard_warp.nii.gz")):
            reg_type = "fnirt"
            interp = interp4fnirt(self.interp)
        else:
            reg_type = "flirt"
            interp = interp4fnirt(self.interp)
        
        infiles = glob(op.join(self.outputdir, "stats", "coef??_*.nii.gz"))
        infiles = glob(op.join(self.outputdir, "stats", "tstat??_*.nii.gz"))
        
        # warp each of the Coef and Tstats
        if reg_type == "flirt":
            for t in ["coef", "tstat"]:
                infnames = glob(op.join(self.outputdir, "stats", "%s??_*.nii.gz" % t))
                for infname in infnames:
                    outprefix = op.basename(infname)
                    cmd = "flirt -in %s -ref reg/standard -applyxfm -init reg/func2standard.mat -interp %s -out reg_standard/%s" % (infname, interp, outprefix)
                    if self.dry_run:
                        self.slog.drycommand(cmd)
                    else:
                        p = self.slog.command(cmd, cwd=self.outputdir)
                        if p.retcode != 0:
                            raise DeconError("problem executing '%s'" % cmd)
        elif reg_type == "fnirt":
            for t in ["coef", "tstat"]:
                infnames = glob(op.join(self.outputdir, "stats", "%s??_*.nii.gz" % t))
                for infname in infnames:
                    outprefix = op.basename(infname)
                    cmd = "applywarp --in=%s --ref=reg/standard --premat=reg/func2highres.mat --warp=reg/highres2standard_warp --interp=%s --out=reg_standard/%s" % (infname, interp, outprefix)
                    if self.dry_run:
                        self.slog.drycommand(cmd)
                    else:
                        p = self.slog.command(cmd, cwd=self.outputdir)
                        if p.retcode != 0:
                            raise DeconError("problem executing '%s'" % cmd)
        
        return
    

class Analyze(object):
    def __init__(self, subject_list, config_file, interp="lin", cores=1, dry_run=False, verbosity=0, 
                 force=False):
        # cores
        os.putenv("OMP_NUM_THREADS", str(cores))
        # config
        self.opts = yaml.load(config_file)
        # args
        self.subject_list = subject_list
        self.interp = interp
        self.dry_run = dry_run
        self.force = force
        # logger
        # note: there is one root logger self.log
        #       and one subject specific logger self.slog
        if verbosity == 0:
            self.loglevel = zlogger.logging.IMPORTANT
        elif verbosity == 1:
            self.loglevel = zlogger.logging.COMMAND_INFO
        else:
            self.loglevel = zlogger.logging.DEBUG
        self.log = zlogger.getLogger("glm_root", level=self.loglevel)
        # setup inputs/outputs and check basics
        self._setup_and_check()
    
    def run(self):
        for s in self.subject_list:
            self.log.title("Running subject: %s" % s)
            errors = []
            
            try:
                # setup variables
                self.log.info("...setting up variables")
                iovars = deepcopy(self.iovars)
                iovars['subject'] = s
                svars = {'subject': s}
                for k in self.iokeys["inputs"]:
                    svars[k] =  zsub_input(iovars[k], iovars)
                for k in self.iokeys["outputs"]:
                    svars[k] =  zsub_output(iovars[k], iovars, self.force)
                for k in self.iokeys["extra"]:
                    svars[k] =  zsub_basic(iovars[k], iovars)
                
                # do it
                analyze = AnalyzeSubject(svars, self.decon, self.reml, self.loglevel, self.interp, self.dry_run)
                analyze.run()
            except IOError as (errno, strerror):
                errors.append(strerror)
            except DeconError as strerror:
                errors.append(strerror)
            except DeconWarning as strerror:
                errors.append(strerror)
            
            if len(errors) > 0:
                self.log.error('You had the following IO errors (skipping to the next subject):')
                for err in errors:
                    self.log.error('\t%s' % err)
            else:
                self.log.debug('You had no errors for this subject!')
        
        return
    
    def _setup_and_check(self):
        self.log.debug("Reading inputs/outputs from config file")
        
        # Basic Check (must have these in the config file)
        self.log.debug("...checking config file")
        check_req(self.opts, ["inputs", "outputdir", "create_model", "run_model"])
        
        # Inputs
        self.log.debug("...getting inputs")
        inputs = self.opts.pop("inputs")
        ## extra stuff
        iovars = inputs.pop("extra", {})
        req_extras = iovars.keys()
        if not isinstance(iovars, dict):
            raise Exception("todo")
        ## funcs, mask, regdir
        req_inputs = ["func", "mask", "regdir"]
        check_req(inputs, req_inputs)
        for k in req_inputs:
            iovars[k] = inputs.pop(k)
        ## check
        for k,v in iovars.iteritems():
            if not isinstance(v, str):
                raise Exception("todo")
        
        # Output
        self.log.debug("...getting outputs")
        outputdir = self.opts.pop("outputdir")
        req_outputs = ["outputdir"]
        if not isinstance(outputdir, str):
            raise Exception("todo")
        iovars["outputdir"] = "%s.reml" % (op.splitext(outputdir)[0])
        
        self.iokeys = {
            "inputs": req_inputs,
            "extra": req_extras,
            "outputs": req_outputs
        }
        
        # 3dDeconvolve
        self.decon = self.opts.pop("create_model")
        
        # 3dREMLfit
        self.reml = self.opts.pop("run_model")
        
        # Save and return
        self.iovars = iovars
        return
    


#####
# Needed Functions
#####

def check_req(d, l):
    """Check that each element of l is in d"""
    for e in l:
        if e not in d:
            raise Exception("todo")

def check_ex(d, l):
    """Check that each element of l is NOT in d"""
    for e in l:
        if e in d:
            raise Exception("todo")

def zsub_basic(template, context):
    """
    Uses string.Template to substitute variables in 'template' with 'context'.
    Does the substitution twice.
    
    template => string
    context => dictionary
    """
    return Template(Template(template).substitute(context)).substitute(context)

def zsub_input(template, context):
    """
    1. Calls zsub_basic
    2. Does glob path expansion
    3. Checks that each resulting path exists
    """
    path = zsub_basic(template, context)
    gpaths = glob(path)
    if len(gpaths) == 0:
        raise IOError("1", "couldn't find input related to '%s'" % path)
    paths = [ op.abspath(op.expanduser(x)) for x in gpaths ]
    if len(paths) == 1:
        return paths[0]
    else:
        return paths


def zsub_output(template, context, force=False):
    """
    1. Calls zsub_basic
    2. Checks the resulting path does not exist but the base directory does exist
    
    TODO: Process the force command here
    """
    path = op.abspath(op.expanduser(zsub_basic(template, context)))
    dpath = op.dirname(path)
    if not force:
        if op.isdir(path) or op.isfile(path):
            raise IOError(2, "output '%s' already exists" % path)
        elif not op.isdir(dpath):
            raise IOError(2, "base directory of output '%s' does not exist" % dpath)
        elif not os.access(dpath, os.W_OK):
            raise IOError(2, "base directory of output '%s' is not writable...fix permissions" % dpath)
    return path

def check_label(label):
    """Display warning if label is greater than 26 characters"""
    l = len(label)
    if l > 24:
        raise DeconWarning(("Label '%s' is too long. " % label) + 
            ("It should have a length of 24 characters or less but has a length of %i." % l))

def interp4fnirt(interp):
    if interp == "lin":
        interp = "trilinear"
    return interp

def interp4flirt(interp):
    if interp == "spline":
        interp = "sinc"
    elif interp == "lin":
        interp = "trilinear"
    elif interp == "nn":
        interp = "nearestneighbour"
    return interp



#####
# Let's Get It On
#####


if __name__ == "__main__":
    main(sys.argv[1:])



# Create 3dDeconvolve here
# todo: create 3dREMLfit (if matrix is specified for 3dDeconvolve)
# todo: using the glt_labels, one should be able to extract the appropriate coef and tstat using 3dcalc
# todo: then register those to standard space using applywarp or flirt command

# todo: later will want to have easy way to apply this to 3dMEMA


# TODO: check for any duplicates    



#####
# Copy data and applywarp?
#####

# loop through glts to get these
#3dcalc -a stats/glts["%s#0_Coef or Tstat" % label] -prefix "coef%02i_%s.nii.gz" % (i, label, ext)

# apply_warp

# also apply easy_thresh?

# remove the original tstats in subject space? or remove the glts?

# in reg_standard also copy over the standard, premat, and fnirt

# can use 3dbucket and 3drefit in order to create a new bucket




#####
# 3dDeconvolve
#####

