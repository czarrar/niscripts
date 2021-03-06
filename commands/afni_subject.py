#!/usr/bin/env python

import argparse, os, sys
sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))

from usage import NiArgumentParser, store_filename, store_input, append_var
from analysis import fromYamlSubject
from zlogger import (LoggerError, LoggerCritical)


#####
# Process user arguments
#####

def create_parser():
    parser = NiArgumentParser(fromfile_prefix_chars='@', argument_default=argparse.SUPPRESS, 
                description="AFNI subject level statistical analysis")
    parser._add_outputs = False
    
    group = parser.add_argument_group('Required')
    group.add_argument('-s', '--subjects', nargs="+", required=True)
    group.add_argument('-c', '--config', action=store_input, check_file=True, required=True, metavar="FILE")
    
    group = parser.add_argument_group('Optional')
    group.add_argument('--var', action="append", type=append_var, dest="vars")
    group.add_argument("--run", nargs="+", dest="run_keys")
    group.add_argument("--decon", action="append_const", const="decon", dest="run_keys")
    group.add_argument("--reml", action="append_const", const="reml", dest="run_keys")
    group.add_argument("--beta-series", action="append_const", const="beta_series", 
                        dest="run_keys")
    group.add_argument("--register-beta-series", action="append_const", 
                        const="register_beta_series", dest="run_keys")
    group.add_argument("--correlate", action="append_const", const="correlate", dest="run_keys")
    group.add_argument("--verbose", action="store_const", const=1, dest="verbosity", default=0)
    group.add_argument("--debug", action="store_const", const=2, dest="verbosity", default=0)
    group.add_argument("--dry-run", action="store_true", default=False)
    group.add_argument("--log", type=store_filename, default=None, metavar="FILE")
    
    return parser

def test_wf(dry_run=True):
    if dry_run:
        arglist = "-c %s/tests/afni_subject.yaml --debug --dry-run -s tb3417" % os.getenv("NISCRIPTS")
    else:
        arglist = "-c %s/tests/afni_subject.yaml--debug -s tb3417" % os.getenv("NISCRIPTS")
    return main(arglist.split())

def main(arglist):
    # Parse
    parser = create_parser()
    args = parser.parse_args(arglist)
    kwargs = vars(args)
    if 'run_keys' not in kwargs:
        kwargs['run_keys'] = ["decon", "reml", "beta_series"]
    if 'vars' not in kwargs:
        template_vars = {}
    else:
        template_vars = dict(kwargs.pop("vars"))
    subjects = kwargs.pop("subjects")
    for subject in subjects:
        try:
            template_vars['subject'] = subject
            fromYamlSubject(user_template_vars=template_vars, **kwargs)
            del template_vars['subject']
        except (LoggerError, LoggerCritical) as err:
            pass


if __name__ == "__main__":
    main(sys.argv[1:])
