#!/usr/bin/env python

import argparse, os, sys
import os.path as op
from glob import glob
from collections import OrderedDict

sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))
import usage
from usage import NiArgumentParser
from report import PreprocReporter
from zlogger import (LoggerError, LoggerCritical)

class store_special(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        labvalues = [ tuple(x.split(":")) for x in values ]
        for x in labvalues:
            if len(x) == 1:
                parser.error("Unrecognized %s. Each argument for '%s' must be a 'label:path'." 
                                % (x, option_string))
        setattr(namespace, self.dest, dict(labvalues))


def create_parser():
    parser = NiArgumentParser(fromfile_prefix_chars="@", 
                              argument_default=argparse.SUPPRESS,
                              description="""Create motion report""")
    parser._add_inputs = False
    parser._add_outputs = False
    
    group = parser.add_argument_group('Required "Options"')
    group.add_argument('-s', '--subjects', nargs="+", required=True)
    group.add_argument('-a', '--anatdirs', nargs="+", action=store_special, required=True)
    group.add_argument('-r', '--regdirs', nargs="+", action=store_special, default=None)
    group.add_argument('-f', '--funcdirs', nargs="+", action=store_special, default=None)
    group.add_argument('-e', '--segdirs', nargs="+", action=store_special, default=None)
    group.add_argument('-n', '--nuidirs', nargs="+", action=store_special, default=None)
    group.add_argument('-o', '--outdir', required=True)
    
    group = parser.add_argument_group('Optional "Options"')
    group.add_argument("--overwrite", action="store_true", default=False)
    group.add_argument("--verbose", action="store_const", const=1, dest="verbosity", default=0)
    group.add_argument("--debug", action="store_const", const=2, dest="verbosity", default=0)
    
    return parser

def main(arglist):
    parser = create_parser()
    args = parser.parse_args(arglist)
    try:
        reporter = PreprocReporter(args.verbosity, args.outdir, overwrite=args.overwrite)
        reporter.setData(args.subjects, args.anatdirs, args.funcdirs, args.regdirs, 
                         args.segdirs, args.nuidirs)
        reporter.run()
    except (LoggerError, LoggerCritical):
        parser.error("Quiting")
    return

if __name__ == "__main__":
    main(sys.argv[1:])
