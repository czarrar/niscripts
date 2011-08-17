class CombineFuncs(SubjectBase):
    """
    Combine Different Functional Runs
    
    Examples
    --------
    >>> from feat_subject import CombineFuncs
    >>> 
    >>> template_vars = {'base': "/Users/zarrar/Projects/tnetworks/output/${subject}/wm", 
                            'subject': 'tb3417'}
    >>> cf = CombineFuncs(2, template_vars)
    >>> 
    >>> infiles = "${base}/run_[1-6]/func_preproc.nii.gz"
    >>> outfunc = "${base}/all_runs/func_preproc.nii.gz"
    >>> cf.setData(infiles, outfunc, mkdir="${base}/all_runs", overwrite=False)
    >>> 
    >>> outmat = "${base}/all_runs/decon_prefeat_mat.1D"
    >>> cf.setDecon(polort="A", censortr="2:37,3:48", outmat=outmat)
    >>> 
    >>> #indirs = "/Users/zarrar/Dropbox/Research/tasks/distwm/stim_files/fsl/tb3417/run_[1-6]"
    >>> #cf.setEvs(indirs=indirs, outdir="${evdir}")
    >>> 
    """
    _logname = "combine_funcs"
    
    #def _combineEVs(self):
    #    if not hasattr(self, "infiles") or len(self.infiles) < 2:
    #        self.log.error("Problem combining EVs...missing or incorrect input files")
    #    many_nvols = [ nibabel.load(f).shape[3] for f in self.infiles ]
    #    n = len(self.infiles)
    #    if n != len(self.evs_outfiles):
    #        self.log.error("Number of functional inputs (%i) does not match number of EV" % n 
    #                        + " inputs (%i)" % len(self.evs_outfiles))
    #    ## new evs
    #    new_nvols = sum(many_nvols)
    #    new_evs = [ np.zeros(new_nvols, self.many_evs[0][i].shape[1]) for i in xrange(n) ]
    #    start_rinds = [0, many_nvols[0]] + [ sum(many_nvols[0:i]) for i in xrange(2,n) ]
    #    end_rinds = (np.array(new_start_rinds) + np.array(many_nvols) - 1).tolist()
    #    ## combine!
    #    for i in xrange(n):
    #        for j,evs in enumerate(self.evs_many):
    #            ev = evs[i]
    #            nvols = many_nvols[i]
    #            if ev.shape[1] == 1:
    #                if ev.shape[0] != nvols:
    #                    self.log.error("# of rows in run %i EV (1-column format) don't match" % i 
    #                                    + " # of volumes for that run (%i)" % nvols)
    #            elif ev.shape[1] == 3:
    #                if max(ev[:,0]) > nvols:
    #                    self.log.error("Start time in run %i EV (3-column format) is longer " % i
    #                                    + "then the length of the run (%i)" % nvols)
    #                ev[:,0] = ev[:,0] + start_rinds[i]                
    #            new_evs[i][start_rinds[i]:end_rinds[i],:] = ev
    #    
    #    self.evs_new = new_evs
    #    return new_evs
    
    def run(self):
        # 3dDeconvolve
        if self.decon_opts:
            self.decon_opts.insert(0, "-input %s" % " ".join(self.infiles))
            self.decon_opts.insert(0, "3dDeconvolve")
            self.log.command(" ".join(self.decon_opts), cwd=op.dirname(self.decon_outmat))
            x = np.loadtxt(self.decon_outmat)
            np.savetxt(self.decon_outmat, x, fmt="%.6f")
        
        ## Combine EVs
        #self._combineEVs()
        #for i in xrange(self.nevs):
        #    np.savetxt(self.evs_outfiles[i], self.evs_new[i])
        
        # Combine Functionals
        cmd = "fslmerge -t %s %s" % (self.outfunc, " ".join(self.infiles))
        self.log.command(cmd, cwd=op.dirname(self.outfunc))
        
        return
    
    def _getInputs(self, infiles, itype):
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
            raise Exception("todo")
        return new_infiles
    
    def setData(self, infiles, outfunc, mkdir=None, overwrite=False):
        """Set inputs and outputs of functional data"""
        self.log.debug("Setting Data")
        
        # Set infiles
        self.infiles = self._getInputs(infiles, 'file')
        
        # Make a directory
        mkdir = self._substitute(mkdir)
        if mkdir and not op.isdir(mkdir):
            self.log.info("Making directory: %s" % mkdir)
            os.mkdir(mkdir)
        
        # Set outfunc
        outfunc = self._substitute(outfunc)
        if op.isfile(outfunc):
            if overwrite:
                self.log.warning("removing output '%s'" % outfunc)
                os.remove(outfunc)
            else:
                raise FeatError("Output '%s' already exists...not overwriting" % outfunc)
        self.outfunc = outfunc
        
        return
    
    def setDecon(self, polort=None, censortr=None, outmat="", tr=None):
        """docstring for setDecon"""
        self.log.debug("Setting 3dDeconvolve")
        decon_opts = []     # add -input later
        
        if polort:
            decon_opts.append("-polort %s" % polort)
        if censortr:
            decon_opts.append("-CENSORTR %s" % censortr)
        if decon_opts:
            if tr:
                decon_opts.append("-force_TR %f" % float(tr))
            if outmat:
                outmat = self._substitute(outmat)
                decon_opts.append("-x1D %s" % outmat)
                self.decon_outmat = outmat
            else:
                raise FeatError("You need to specify the outmat option to run 3dDeconvolve")
            decon_opts.append("-x1D_stop")
        
        self.decon_opts = decon_opts
        return
    
    #def setEvs(self, indirs, outdir):
    #    """docstring for setEvs"""
    #    self.log.debug("Setting EVs stuff")
    #    
    #    # Inputs
    #    self.log.debug("...inputs")
    #    indirs = self._getInputs(indirs, 'directory')
    #    many_infiles = [ [ op.join(x, y) for y in os.listdir(x) ] for x in indirs ]
    #    many_evs = [ [ np.loadtxt(f) for f in fs ] for fs in many_infiles ]
    #    many_ncols = [ [ ev.shape[1] for ev in evs ] for evs in many_evs]
    #    ## check
    #    ref_infiles = many_infiles[0]
    #    ref_ncols = many_ncols[0]
    #    for i,dirname in enumerate(indirs):
    #        infiles = many_infiles[i]
    #        ncols = many_ncols[i]
    #        if ref_infiles != infiles:
    #            raise FeatError("Files in EV dir '%s' do not match those in '%s'" % (ref_indir, 
    #                                                                                    dirname))
    #        for j,ncol in enumerate(ncols):
    #            if ref_ncols[j] != ncol:
    #                raise FeatError("Number of columns in '%s' do not match those in '%s'" % 
    #                                    (ref_infiles[j], infiles[j]))
    #            elif ncol != 1 and ncol !=3:
    #                raise FeatError("EV file '%s' must be 1 or 3 column format" % ref_infiles[j])
    #    ## save
    #    self.evs_many = many_evs
    #    
    #    # Output
    #    self.log.debug("...outputs")
    #    outdir = self._substitute(outdir)
    #    if not op.isdir(outdir):
    #        self.log.info("Creating EV output directory" % outdir)
    #        os.mkdir(outdir)
    #    self.evs_outfiles = [ op.join(outdir, op.basename(p)) for p in ref_infiles ]
    #    self.nevs = len(self.evs_outfiles)
    #    
    #    return

