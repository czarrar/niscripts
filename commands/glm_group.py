#!/usr/bin/env python

# glob(/path/to/subject.reml)
# go through each reml directory and find the coefs and tstats
# find ones that overlap, give warning for any that don't

class AnalyzeGroup(object):
    def __init__(self, dry_run=False, verbosity=0, force=False):
        self.dry_run = dry_run
        self.force = force
        if verbosity == 0:
            self.loglevel = zlogger.logging.IMPORTANT
        elif verbosity == 1:
            self.loglevel = zlogger.logging.COMMAND_INFO
        else:
            self.loglevel = zlogger.logging.DEBUG
        self.rlog = zlogger.getLogger("glm_root", level=self.loglevel)
    
    
    def _setupMema(self, sets, model_type, **opts):
        """
        Creates command-line stuff for 3dMEMA
        
        prefix: output prefix
        sets: a dict with each set of groups, values of dict must be list of tuples with len=3
        model_type: can be 'one-sample', 'paired', or 'two-sample'
        **opts: all the other opt=args, note for opts without args, use True/False for arg
        """
        self.rlog.debug("Setup Mema")
        
        mema_opts = []
        
        # Prefix
        mema_opts.append("-prefix stat")
        
        # Check/set groups
        groups = sets.keys()
        ngroups = len(sets)
        if ngroups == 1 and model_type != "one-sample":
            raise MemaError("One group set given...model type must be 'one-sample'")
        if ngroups == 2 and model_type == "one-sample":
            raise MemaError("Two group sets given...model type must be 'paired' or 'two-sample'")
        if ngroups > 2:
            raise MemaError("Can only have a maximum of two group sets")
        
        # Setup conditions/groups
        if ngroups == 2:
            if model_type == "paired":
                mema_opts.append("-conditions %s", " ".join(groups))
            elif model_type == "two-sample":
                mema_opts.append("-groups %s", " ".join(groups))
        
        # Setup sets
        set_opts = []
        for k,vs in sets:
            vs = [ " ".join(x) for x in vs if len(x) == 3 ]
            mema_opts.append("-set %s %s" % (k, " ".join(vs)))
        
        # TRUE/FALSE ones
        boolean_no_opts = ["HKtest", "model_outliers", "residual_Z"]
        boolean_un_opts = ["equal_variance"]
        boolean_opts = boolean_no_opts + boolean_un_opts
        for k in boolean_opts:
            if k in opts:
                v = opts.pop(k)
                if not isinstance(v, bool):
                    raise MemaError("Option '%s' must have argument being True or False" % k)
                elif v:
                    mema_opts.append("-%s" % k)
                elif k in boolean_no_opts:
                    mema_opts.append("-no_%s" % k)
                elif k in boolean_un_opts:
                    mema_opts.append("-un%s" % k)
        
        # Everything else
        for k,v in opts:
            mema_opts.append("-%s %s" % (k,v))
        
        mema_opts.insert(0, "3dMEMA")
        self.mema_opts = mema_opts
        
        return self.mema_opts
    
    def _setupManySets(self, inputdir, subjects, contrasts):
        """
        Create the proper input for sets to go into setupMema
        
        inputdir: a string with ${subject} variable to represent each input subject reml dir
        subjects: list of subjects
        set_names: list of group names
        contrasts: list of different contrasts present for each subject
        
        Returns tuple(names, sets)
        """
        self.rlog.debug("Setup Many Sets")
        
        if inputdir.keys() != contrasts.keys():
            raise MemaError("inputdir and contrasts must have same element names")
        names = inputdir.keys()
        
        sets = {}
        for c in contrasts:
            self.rlog.debug("...contrast: %s" % c)
            sets[c] = {}
            for n in names:
                self.rlog.debug("......group: %s" % n)
                tmp = []
                for s in subjects:
                    self.rlog.debug(".........subject: %s" % s)
                    sdir = Template(inputdir[n]).context(subject=s)
                    if not op.isdir(sdir):
                        raise MemaError("Directory '%s' does not exist" % sdir)
                    gpath1 = op.join(sdir, "reg_standard", "coefs+tlrc'[%s]'" % con)
                    gpath2 = op.join(sdir, "reg_standard", "tstats+tlrc'[%s]'" % con)
                    tmp.append((s, gpath1, gpath2))
                sets[c][n] = tmp
        
        return sets
    
    def setupStats(self, inputdir, subjects, contrasts, model_type, **opts):
        self.rlog.info("Setting up stats commands")
        
        check_ex(opts, ["mask", "prefix"])
        
        ncontrasts = len(contrasts)
        self.contrasts = contrasts
        
        many_sets = self._setupManySets(inputdir, subjects, contrasts)
        self.many_mema_opts = [ 
            self._setupMema(many_sets[cname], model_type, **opts) for cname in contrasts 
        ]
        
        return self.many_mema_opts
    
    # todo: add some intermediate thing that will copy the stats and create appropriate stuff
    
    def setupPostStats(self):
        self.rlog.info("Setting up post-stats commands")
        
        cope%02i_*/zstat1 mask self.voxthresh self.clusthresh standard grot
        
        "easythresh stats/zstat1 mask 2.3 0.01 example_func grot"
    
    def fromYaml(self, config_file, subjects):
        self.rootlog.debug("Loading info from yaml file")
        opts = yaml.load(config_file)
        check_req(opts, ["inputdirs", "outputdir", "stats", "post-stats"])
        
        # Input/Output
        inputdir = opts.pop("inputdir")
        outputdir = opts.pop("outputdir")
        if op.isdir(outputdir):
            raise MemaError("Output directory '%s' already exists" % outputdir)
        
        # Stats
        stats = opts.pop("stats")
        self.setupStats(stats)
        
        # Post-Stats
        post_stats = opts.pop("post-stats")
        self.setupPostStats(post_stats)
        
        # call setup
        self.setupCommand(inputdir, outputdir, **stats)
        
        # 
        
        
        ## Set
        self.log.debug("...getting set names")
        sets = self.opts.pop("sets")
        if isinstance(sets, str):
            ngroups = 1
        elif isinstance(sets, list):
            if len(sets) == 1:
                sets = sets[0]
                ngroups = 1
            elif len(sets) == 2:
                ngroups = 2
            else:
                raise MemaError("can only have up to 2 sets")
        self.mema['set'] = sets
        self.mema['ngroups'] = ngroups
        
        ## Contrasts
        self.log.debug("...getting contrasts")
        contrasts = self.opts.pop("contrasts")
        if not isinstance(contrasts, list):
            raise MemaError("contrasts must be a list")
        self.mema['contrasts'] = contrasts
        
        
        if t.find("${subject}") == -1 or t.find("$subject") == -1:
            raise MemaError("inputdirs must contain ${subject} or $subject variable")
        self.inputdirs = [ zsub_input(t, {'subject': s}) for s in subject_list ]
        
        t = "%s.mema" % (op.splitext(t)[0])
        self.outputdir = zsub_output(t, {})
        
    

    
    
    
    
    def _setup_and_check(self, subject_list):
        self.log.info("Getting stuff from config file")
        
        error = None
        try:
            # Basic Check (must have these in the config file)
            self.log.debug("...checking config file")
            check_req(self.opts, ["sets", "contrasts", "inputdirs", "outputdir"])
            check_ex(self.opts, ["mask", "prefix"])
        
            # Inputs
            self.log.debug("...getting inputs")
            t = self.opts.pop("inputdirs")
            if not isinstance(t, str):
                raise MemaError("inputdirs must be a string")
            if t.find("${subject}") == -1 or t.find("$subject") == -1:
                raise MemaError("inputdirs must contain ${subject} or $subject variable")
            self.inputdirs = [ zsub_input(t, {'subject': s}) for s in subject_list ]
        
            # Output
            self.log.debug("...getting output")
            t = self.opts.pop("outputdir")
            if not isinstance(t, str):
                raise MemaError("outputdir must be a string")
            t = "%s.mema" % (op.splitext(t)[0])
            self.outputdir = zsub_output(t, {})
        
            # MEMA options
            self.mema = {}
            
            ## Set
            self.log.debug("...getting set names")
            sets = self.opts.pop("sets")
            if isinstance(sets, str):
                ngroups = 1
            elif isinstance(sets, list):
                if len(sets) == 1:
                    sets = sets[0]
                    ngroups = 1
                elif len(sets) == 2:
                    ngroups = 2
                else:
                    raise MemaError("can only have up to 2 sets")
            self.mema['set'] = sets
            self.mema['ngroups'] = ngroups
            
            ## Contrasts
            self.log.debug("...getting contrasts")
            contrasts = self.opts.pop("contrasts")
            if not isinstance(contrasts, list):
                raise MemaError("contrasts must be a list")
            self.mema['contrasts'] = contrasts
            
            ## other stuff
            self.mema.update(self.opts)
        except IOError as (errno, strerror):
            error = strerror
        except DeconError as strerror:
            error = strerror
        
        if error:
            self.log.error('Error detected: %s' % error)
            raise SystemExit(1)
        
        return
    
    def _create_output(self):
        os.mkdir(self.outputdir)
    
    def run(self):
        self.log.title("Running MEMA")
        
        # Gather inputs
        
        # Create output
        self.log.info("creating output")
        self._create_output()
        


#####
# Needed Functions
#####

def check_req(d, l):
    """Check that each element of l is in d"""
    for e in l:
        if e not in d:
            raise MemaError("%s must be in options" % e)

def check_ex(d, l):
    """Check that each element of l is NOT in d"""
    for e in l:
        if e in d:
            raise MemaError("%s must not be in options" % e)

def zglob_one(path):
    gpath = glob(path)
    if len(gpath) == 0:
        raise MemaError("Nothing found for %s" % path)
    elif len(gpath) > 1:
        raise MemaError("%i things found for %s but should have been 1" % (len(gpath), path))
    return gpath[0]

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
