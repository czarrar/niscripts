#!/usr/bin/env python

import argparse, os, sys
sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))

from usage import NiArgumentParser, store_filename, store_input
from analysis import fromYamlSubject
from zlogger import (LoggerError, LoggerCritical)

#####
# Process user arguments
#####

def create_parser():
    parser = NiArgumentParser(fromfile_prefix_chars='@', argument_default=argparse.SUPPRESS, 
                description="FSL subject level statistical analysis (creates model)")
    parser._add_outputs = False
    
    group = parser.add_argument_group('Required')
    group.add_argument('-s', '--subjects', nargs="+", required=True)
    group.add_argument('-c', '--config', action=store_input, check_file=True, required=True)
    
    group = parser.add_argument_group('Optional')
    group.add_argument("--combine", action="append_const", const="combine", dest="run_keys")
    group.add_argument("--res-decon", action="append_const", const="res_decon", 
                        dest="run_keys")
    group.add_argument("--fsf", action="append_const", const="fsf", dest="run_keys")
    group.add_argument("--feat", action="append_const", const="feat", dest="run_keys")
    group.add_argument("--regress", action="append_const", const="regress", dest="run_keys")
    group.add_argument("--verbose", action="store_const", const=1, dest="verbosity", default=0)
    group.add_argument("--debug", action="store_const", const=2, dest="verbosity", default=0)
    group.add_argument("--dry-run", action="store_true", default=False)
    group.add_argument("--log", type=store_filename, default=None, metavar="FILE")
    
    return parser

def test_wf(dry_run=True):
    if dry_run:
        arglist = "-c %s/tests/feat_subject.yaml --dry-run --debug -s tb3417" % os.getenv("NISCRIPTS")
    else:
        arglist = "-c %s/tests/feat_subject.yaml--debug -s tb3417" % os.getenv("NISCRIPTS")
    return main(arglist.split())

def main(arglist):
    # Parse
    parser = create_parser()
    args = parser.parse_args(arglist)
    kwargs = vars(args)
    if 'run_keys' not in kwargs and not kwargs['run_keys']:
        kwargs['run_keys'] = ["res_decon", "fsf", "feat", "regress"]
    subjects = kwargs.pop("subjects")
    for subject in subjects:
        try:
            fromYamlSubject(subject=subject, **kwargs)
        except (LoggerError, LoggerCritical) as err:
            pass


if __name__ == "__main__":
    main(sys.argv[1:])
