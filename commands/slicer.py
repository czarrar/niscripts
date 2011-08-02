#!/usr/bin/env python

# Author: Zarrar Shehzad
# Version: 0.1
# Date: July 4, 2011

# This is a wrapper around FSL's slicer and overlay programs
# note: several options from these two programs aren't supported here

import os, sys
sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))

from execute import *
from tempfile import mkstemp, mkdtemp
from math import ceil
import re, argparse, string

###
# Process command-line inputs
###

parser = argparse.ArgumentParser(
    description="""
        Generate PNG images from 2D brain slices. Can optionally overlay 
        statistical images or overlay an outline of a brain image. This is a 
        wrapper around FSL's overlay, slicer, and pngappend.
    """)


# ACTIONS USED LATER

class store_input(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if not os.path.isfile(value):
            parser.error("Input file '%s' does not exist" % value)
        setattr(namespace, self.dest, os.path.abspath(value))
        # also store input dimensions
        tmp = sh.fslinfo(value).stdout
        info = dict([ tuple(re.split("[\ ]+", x)) for x in tmp.split("\n") ])
        namespace._dim = {
            'x': int(info["dim1"]), 
            'y': int(info["dim2"]), 
            'z': int(info["dim3"]),
        }
    

class store_filename(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if not os.path.isfile(value):
            parser.error("File '%s' does not exist" % value)
        setattr(namespace, self.dest, os.path.abspath(value))
    

class store_overlay(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not os.path.isfile(values[0]):
            parser.error("Overlay file '%s' does not exist" % values[0])
        values[0] = os.path.abspath(values[0])
        setattr(namespace, self.dest, values)
    

class store_output(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        value = re.sub("[.].{1,4}$", "", value)
        value = value + ".png"
        setattr(namespace, self.dest, os.path.abspath(value))
    

class append_to_main(argparse.Action):
    to_fsl = {
        '--slice-labels': '-L',
        '--lut': '-l',
        '--scale': '-s',
        '--intensity': '-i',
        '--no-lr-labels': '-u'
    }
    
    def __call__(self, parser, namespace, values, option_string):
        if not hasattr(namespace, 'main'): namespace.main = []
        if isinstance(values, list):
            values = string.join(values, " ")
        if values:
            namespace.main.append('%s %s' % (self.to_fsl[option_string], 
                values))
        else:
            namespace.main.append('%s' % self.to_fsl[option_string])
    

# GENERAL INPUT/OUTPUT

group = parser.add_argument_group('Positional Arguments')
group.add_argument('input', action=store_input, help="Input 3D brain image")
group.add_argument('output', action=store_output, help="Output 2D PNG image")


# OVERLAYS

group = parser.add_argument_group('Overlay Inputs')

## FSL slicer
group.add_argument('--edge-overlay', action=store_filename, metavar="FILE", 
    default=None, help="""An outline of a brain image is overlaid on the input
    (useful for checking registration). Note: the overlay options don't apply
    to this input.""")

## FSL overlay
group.add_argument('--overlay', action=store_overlay, nargs=3, default=None,
    metavar=("file", "min", "max"), help="image to overlay on input")
exclusive_group = group.add_mutually_exclusive_group()
exclusive_group.add_argument('--overlay2', action=store_overlay, nargs=3, 
    default=None, help="2nd image to overlay on input", metavar=("file", "min", 
    "max"), )
exclusive_group.add_argument("--show-negative", action="store_true", 
    default=False, help="""Will take the 1st overlay and threshold it to have a 
    negative range (this would be useful if you wanted to display positive & 
    negative thresholded maps). Note: you can't specify both the --overlay2 and
    this option.""")


# OVERLAY OPTIONS

group = parser.add_argument_group('Overlay Options', 
    'Note: these only apply to --overlay and --overlay2')

group.add_argument('-t', '--transparency', action="store_const", default="0",
    const="1", help="make overlay colors semi-transparent")
group.add_argument('--checkerboard', action="store_true", default=False,
    help="use checkerboard mask for overlay")
group.add_argument('--background-range', nargs=2, metavar=("min", "max"),
    default=["-a"], help="setting range of underlay, only used when overlay " + 
    "option has been given (default: automatic estimation)")


# MAIN SLICER OPTIONS

group = parser.add_argument_group('Main Image Options')

## --slice-labels (-L)
group.add_argument('--slice-labels', action=append_to_main, nargs=0, 
    default=argparse.SUPPRESS, dest="main")

## -l (-l)
group.add_argument('--lut', action=append_to_main, default=argparse.SUPPRESS, 
    metavar="color", help="use a different colour map from that specified in" +
    " the header", dest="main")

## -s/--scale (-s)
group.add_argument('--scale', action=append_to_main, default=argparse.SUPPRESS, 
    metavar="positive-#", help="relative size of each slice", dest="main")

## --intensity (-i)
group.add_argument('--intensity', action=append_to_main, nargs=2, 
    default=argparse.SUPPRESS, metavar=("min", "max"), help="specify intensity"+  
    " min and max for display range", dest="main")

## --no-lr-labels (-u)
group.add_argument('--no-lr-labels', action=append_to_main, nargs=0, 
    default=argparse.SUPPRESS, dest="main")


# OUTPUT SLICER OPTIONS

group = parser.add_argument_group('Output Image Options')

## slice type
group.add_argument('-s', '--slice', required=True, metavar="slice_name", 
    help="type of slice; can be: a/axial, c/coronal, or s/sagittal")

## image width
group.add_argument('-w', '--width', required=True, type=int, metavar="slices", 
    help="image width in # of slices")

## image height or slice step
exclusive_group = group.add_mutually_exclusive_group(required=True)
exclusive_group.add_argument('-l', '--height', default=argparse.SUPPRESS, 
    type=int, help="image height in # of slices", metavar="slices")
exclusive_group.add_argument('-e', '--slice-every', dest="sample", type=int,
    default=argparse.SUPPRESS, help="include every X # of slice", 
    metavar="number")


# OTHER STUFF

group = parser.add_argument_group('Miscellaneous Options')

group.add_argument('--force', action='store_true', default=False, 
    help="overwrite output if it exists")
group.add_argument('-v', '--verbose', action="store_true", default=False)
group.add_argument('--dry-run', action="store_true", default=False, 
    help="Won't execute anything but prints what would have been executed")
group.add_argument('--version', action='version', version="%(prog)s 0.1")



# FSL's Overlay
def compile_overlay_args(parser, args):
    # Start
    overlay_args = []
    
    # Transparency?
    overlay_args.append(args.transparency)
    
    # Output: floating point 
    overlay_args.append("0")
    
    # Checkerboard
    if args.checkerboard:
        overlay_args.append("-c")
    
    # Underlay
    overlay_args.append(args.input)
    
    # Underlay Range
    overlay_args.append(string.join(args.background_range, " "))
    
    # Overlay
    overlay_args.append(string.join(args.overlay, " "))
    
    # Overlay 2
    if args.overlay2:
        overlay_args.append(string.join(args.overlay2, " "))
    elif args.show_negative:
        args.overlay2 = args.overlay
        args.overlay2[1] = "-%s" % args.overlay2[1]
        args.overlay2[2] = "-%s" % args.overlay2[2]
        overlay_args.append(string.join(args.overlay2, " "))
    
    # Temporary Output
    ## TODO: have suffix determined from FSLOUTPUTTYPE
    outfname = "overlay.nii.gz"
    overlay_args.append(outfname)
    
    return (string.join(overlay_args, " "), outfname)


def compile_slicer_args(parser, args):
    # Start
    slicer_args = []
    
    # Input
    slicer_args.append(args.input)
    
    # Input2?
    if args.edge_overlay:
        slicer_args.append(args.edge_overlay)
    
    # Main Options
    if hasattr(args, "main"):
        slicer_args.extend(args.main)
    
    # Output Options
    
    ## slice name
    slice_name_lookup = {
        'sagittal': 'x', 's': 'x', 
        'coronal': 'y', 'c': 'y',
        'axial': 'z', 'a': 'z',
    }
    if args.slice in slice_name_lookup:
        skey = slice_name_lookup[args.slice]
    else:
        parser.error("slice name must be one of the following: \n%s" %
            string.join(slice_name_lookup.keys(), ", "))
    
    ## total number of slices
    if args.height:
        nslices = args.width * args.height
    elif args.sample:
        nslices = round(float(args._dim[skey])/args.sample)
    
    ## list of absolute slice #s, etc
    step = int(ceil(float(args._dim[skey])/nslices))
    slice_numbers = range(0, -args._dim[skey], -step)
    nslices = len(slice_numbers)
    
    # save for later
    args.height = int(ceil(float(nslices)/args.width))
    args._nslices = nslices
    
    ## temporary directory
    slice_fnames = [ 
        "slice_%s%i.png" % (skey, i) for i in range(1, nslices+1)
    ]
    
    ## add options
    for n,f in zip(slice_numbers, slice_fnames):
        slicer_args.append('-%s %i %s' % (skey, n, f))
    
    return (string.join(slicer_args, " "), slice_fnames)


def compile_pngappend_args(parser, args, slice_fnames):
    # input images to append
    w = args.width; h = args.height
    pngappend_args = string.join([ 
        string.join([ 
            slice_fnames[i] for i in range(x*w, x*w+w) if i < args._nslices
        ], " + ")
        for x in range(h)
    ], " - ")
    
    # output
    pngappend_args += " " + args.output
    
    return pngappend_args


if __name__ == "__main__":
    args = parser.parse_args()
    
    # check output
    if os.path.isfile(args.output) and not args.force:
        parser.exit("Output file '%s' already exists, use --force if you want to overwrite it" % args.output)
        
    # TODO: check for required programs on path
    
    curdir = os.getcwd()
    tmpdir = None; tmp_input = None; slice_fnames = []
    try:
        # create temporary directory
        tmpdir = mkdtemp(prefix='slicer_tmp')
        if args.verbose or args.dry_run:
            print '\nmkdir %s' % tmpdir
            print 'cd %s' % tmpdir
        os.chdir(tmpdir)
        
        # overlay
        if args.overlay:
            (overlay_args, tmp_input) = compile_overlay_args(parser, args)
            if args.verbose or args.dry_run:
                print "\noverlay %s" % overlay_args
            if not args.dry_run:
                result = sh.overlay(overlay_args)
                if result.retcode:
                    parser.exit(3, "error running overlay: \n%s\n" %
                        result.stderr)
            args.input = tmp_input
        
        # slicer
        (slicer_args, slice_fnames) = compile_slicer_args(parser, args)
        if args.verbose or args.dry_run:
            print "\nslicer %s" % slicer_args
        if not args.dry_run:
            result = sh.slicer(slicer_args)
            if result.retcode:
                parser.exit(3, "error running slicer: \n%s\n" % result.stderr)
        
        # pngappend
        pngappend_args = compile_pngappend_args(parser, args, slice_fnames)
        if args.verbose or args.dry_run:
            print "\npngappend %s" % pngappend_args
        if not args.dry_run:
            result = sh.pngappend(pngappend_args)
            if result.retcode:
                parser.exit(3, "error running pngappend: \n%s\n" % 
                    result.stderr)
    finally:
        if args.verbose:
            print "\n...cleaning up\n"
        # Remove temporary directories/files
        if tmp_input:
            os.remove(tmp_input)
        if not args.dry_run:
            for f in slice_fnames:
                os.remove(f)
        if tmpdir:
            os.rmdir(tmpdir)
        os.chdir(curdir)
    
    sys.exit()
