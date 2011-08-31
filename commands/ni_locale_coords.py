#!/usr/bin/env python

import argparse, os, sys
import os.path as op
import numpy as np

import nipype.interfaces.fsl as fsl # fsl
standard_image = fsl.Info.standard_image

sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))
import zlogger
from usage import NiArgumentParser, store_filename, store_input, store_output, append_var
from zlogger import (LoggerError, LoggerCritical)
from analysis.base import Base, get_loglevel, create_logger
from execute import Process

class CoordinateSubject(Base):
    """docstring for CoordinateSubject"""
    
    _logname = "coordinate_subject"
    _output_suffix = {
        "peaks": "table.txt",
        "filt": "table_filt.1D",
        "coord": "peak.1D",
        "roi": "roi.nii.gz"
    }
    
    def __init__(self, *args, **kwargs):
        super(CoordinateSubject, self).__init__(*args, **kwargs)
        self.log.info("Starting CoordinateSubject")
        self._isset_data = self._isset_options = False
        self.infile = self.inmask = self.thresh = self.refcoord = self.roi_rad = self.orient = None
        self.outputs = {}
        self.cmd_extrema = self.cmd_filt = self.cmd_roi = ""
        self.cmds_filt = []
        return
    
    def compile(self):
        self.log.info("Compiling")
        if not self._isset_data and not self._isset_options:
            self.log.fatal("Must set data and options")
        self.cmd_extrema = "3dExtrema -data_thr %.3f -mask_file %s -volume %s" % (
                                self.thresh, self.inmask, self.infile)
        
        self.cmd_filt = "tail -n +11 %s" % self.outputs['peaks']
        
        self.cmd_roi = "3dUndump -prefix %s -master %s -srad %i -orient %s -xyz %s" % (
                            self.outputs['roi'], self.std, self.roi_rad, self.orient, 
                            self.outputs['coord'])
        
    
    def _run_closest_coord(self):
        if not self._isset_data:
            self.log.fatal("Must set data")
        x = np.loadtxt(self.outputs["filt"])
        c = x[:,2:5]
        ref = self.refcoord
        self.log.debug("...reference (RAI): %s" % " ".join([ str(x) for x in ref ]))
        dists = np.sqrt((c[:,0]-ref[0])**2 + (c[:,1]-ref[1])**2 + (c[:,2]-ref[2])**2)
        closest_coord = [ int(x) for x in c[np.argmin(dists),:] ]
        self.log.debug("...closest peak (RAI): %s" % " ".join([ str(x) for x in closest_coord]))
        closest_coord[0] = closest_coord[0]*-1
        closest_coord[1] = closest_coord[1]*-1
        disp = "%i %i %i" % (closest_coord[0], closest_coord[1], closest_coord[2])
        self.log.important(disp)
        f = file(self.outputs['coord'], 'w')
        f.write(disp + "\n")
        f.close()
        return closest_coord
    
    def run(self):
        self.compile()
        
        self.log.info("Running")
        if not self.cmd_extrema and not self.cmds_filt:
            self.log.fatal("Must compile")
        
        if self.dry_run:
            self.log.drycommand(self.cmd_extrema + " > %s" % self.outputs['peaks'])
            self.log.drycommand("%s | %s > %s" % (self.cmds_filt[0], self.cmds_filt[1], 
                                                  self.outputs['filt']))
            if self.outroi:
                self.log.drycommand(self.cmd_roi)
        else:
            self.log.debug("...extrema")
            f = file(self.outputs["peaks"], 'w')
            self.log.drycommand("%s > %s" % (self.cmd_extrema, self.outputs["peaks"]))
            p = Process(self.cmd_extrema, stdout=f)
            if p.retcode != 0:
                self.log.error("problem with 3dExtrema")
            if p.stderr:
                print p.stderr
            f.close()
            
            self.log.debug("...filter")
            self.log.drycommand("%s > %s" % (self.cmd_filt, self.outputs['filt']))
            f = file(self.outputs["filt"], 'w')
            p = Process(self.cmd_filt, stdout=f)
            if p.retcode != 0:
                self.log.error("Error using head/tail commands")
            if p.stderr:
                print p.stderr
            f.close()
            
            self.log.debug("...closest coordinate")
            res = self._run_closest_coord()
            
            if self.outroi:
                self.log.debug("...save ROI")
                self.log.command(self.cmd_roi)
            
            return res
    
    def setData(self, infile, inmask, std, outprefix, refcoord, outroi=False, overwrite=False):
        self.log.info("Setting up data")
        
        # Inputs
        self.infile = self._check_infile(infile, desc="input file", substitute=True)
        self.inmask = self._check_infile(inmask, desc="input mask", substitute=True)
        self.std = self._check_infile(std, desc="standard image", substitute=True)
        
        # Outputs
        self.outputs = {}
        outprefix = self._substitute(outprefix)
        for k,v in self._output_suffix.iteritems():
            out = "%s_%s" % (outprefix, v)
            out = self._check_outfile(out, desc="%s file" % k, substitute=True, 
                                      overwrite=overwrite)
            self.outputs[k] = out
        if not op.isdir(op.dirname(outprefix)):
            self.log.info("Creating base directory '%s' for output" % op.dirname(outprefix))
            os.mkdir(op.dirname(outprefix))
        self.outroi = outroi
        
        # Coordinate
        if len(refcoord) != 3:
            self.log.fatal("Must have x,y,z for reference coordinate")
        self.refcoord = np.array([ float(c) for c in refcoord ])
        self.refcoord[0] = self.refcoord[0]*-1
        self.refcoord[1] = self.refcoord[1]*-1
        
        self._isset_data = True
        return
    
    def setOptions(self, thresh=2.3, roi_rad=2, orient="LPI"):
        self.log.info("Setting options")
        
        self.thresh = 2.3
        self.roi_rad = roi_rad
        self.orient = orient
        
        self._isset_options = True
        return
    

def create_parser():
    parser = NiArgumentParser(fromfile_prefix_chars='@', argument_default=argparse.SUPPRESS, 
                description="Find peak coordinate for each subject that is closest to a given" \
                            " reference coordinate")
    parser._add_inputs = False
    parser._add_outputs = False
    
    group = parser.add_argument_group('Required')
    group.add_argument('-i', '--input', required=True, metavar="FILE")
    group.add_argument('-m', '--mask', required=True, metavar="FILE")
    group.add_argument('-o', '--prefix', required=True, metavar="FILE PREFIX")
    group.add_argument('-c', '--coord', nargs=3, type=int, required=True, metavar=('x','y','z'))
    
    group = parser.add_argument_group('Command Options')
    group.add_argument("-t", "--thresh", type=float, default=2.3, metavar="THRESHOLD", help="default: %(default)s")
    group.add_argument('--std', default=standard_image("MNI152_T1_2mm.nii.gz"), metavar="FILE", help="default: %(default)s")
    group.add_argument("--orient", default="LPI", metavar="XXX", help="default: %(default)s")
    group.add_argument('-r', "--radius", type=int, default=1, metavar="mm", help="default: %(default)s")
    
    group = parser.add_argument_group('I/O Options')
    group.add_argument('--roi', action="store_true", default=False, help="create ROI")
    group.add_argument('-s', '--subjects', nargs="+", metavar="ID")
    group.add_argument('--var', action="append", type=append_var, dest="vars")
    group.add_argument('--overwrite', action="store_true", default=False, help="default: %(default)s")
    group.add_argument("--verbose", action="store_const", const=1, dest="verbosity", default=0)
    group.add_argument("--debug", action="store_const", const=2, dest="verbosity", default=0)
    group.add_argument("--dry-run", action="store_true", default=False, help="default: %(default)s")
    group.add_argument("--log", type=store_filename, default=None, metavar="FILE", help="default: %(default)s")
    
    return parser

def run_subject(log, args, template_vars):
    ## setup
    cs = CoordinateSubject(args.verbosity, deepcopy(template_vars), args.dry_run, logger=log)
    cs.setData(args.input, args.mask, args.std, args.prefix, args.coord, args.roi, args.overwrite)
    cs.setOptions(args.thresh, args.radius, args.orient)
    ## run
    cs.run()
    del cs
    return

def main(arglist):
    # Parse
    parser = create_parser()
    args = parser.parse_args(arglist)
    kwargs = vars(args)
    
    # Template vars
    if 'vars' not in kwargs:
        template_vars = {}
    else:
        template_vars = dict(kwargs.pop("vars"))
    
    # Logger
    loglevel = get_loglevel(args.verbosity)
    log = create_logger("coordinate_subjects", loglevel, args.log)
    
    # Loop through subjects
    if 'subjects' in kwargs:
        subjects = kwargs.pop("subjects")
        for subject in subjects:
            log.title("Subject: %s" % subject)
            template_vars['subject'] = subject
            run_subject(log, args, template_vars)
            del template_vars['subject']
    else:
        run_subject(log, args, template_vars)
    return

if __name__ == "__main__":
    main(sys.argv[1:])
