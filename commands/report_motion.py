#!/usr/bin/env python

import argparse, os, sys
import os.path as op
from glob import glob

sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))
import usage
from usage import NiArgumentParser
from report import MotionReporter

def create_parser():
    parser = NiArgumentParser(fromfile_prefix_chars="@", 
                              argument_default=argparse.SUPPRESS,
                              description="""Create motion report""")
    parser._add_inputs = False
    parser._add_outputs = False
    
    group = parser.add_argument_group('Required "Options"')
    group.add_argument('-b', '--base-dir', action=usage.store_directory, required=True)
    group.add_argument('-s', '--subjects', nargs="+", required=True)
    group.add_argument('-r', '--run-dirs', required=True)
    group.add_argument('-o', '--output', required=True)
    
    group = parser.add_argument_group('Optional "Options"')
    group.add_argument("--verbose", action="store_const", const=1, dest="verbosity", default=0)
    group.add_argument("--debug", action="store_const", const=2, dest="verbosity", default=0)
    
    return parser

def main(arglist):
    parser = create_parser()
    args = parser.parse_args(arglist)
    sinfo = dict([ (s, op.join(args.base_dir, s, args.run_dirs)) for s in args.subjects ])
    
    reporter = MotionReporter(args.verbosity)
    reporter.setData(sinfo, args.output)
    reporter.run()
    
    return

if __name__ == "__main__":
    main(sys.argv[1:])
