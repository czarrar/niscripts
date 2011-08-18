#!/usr/bin/env python

# TODO: create  

import argparse, os, sys
sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))

import usage
from analysis import fromYamlSubject, FeatSubject
from zlogger import (LoggerError, LoggerCritical)

def feat(name, workingdir, subject_list, **kwards):
    
    wf = pe.Workflow(name=name)
    wf.base_dir = workingdir
    
    subinfo = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name='subinfo', 
                        iterables=('subject_id', subject_list))
    
    featsubject = FeatSubject()
    for k in kwards:
        setattr(featsubject.inputs, k, v)
    fsnode = pe.Node(interface=featsubject, name="featpy")
    
    wf.connect(subinfo, 'subject_id', fsnode, 'subject')
    
    return wf


#####
# Process user arguments
#####

class FeatParser(usage.NiParser):
    def _create_parser(self, *args, **kwrds):
        """Create command-line interface"""
        parser = usage.NiArgumentParser(fromfile_prefix_chars="@", 
                                  argument_default=argparse.SUPPRESS,
                                  description="""Wrapper around feat_subject_worker.py that parallelizes processing of different subjects.""")
        parser._add_inputs = False
        parser._add_outputs = False
        
        group = parser.add_argument_group('General Options')
        group.add_argument('-s', '--subjects', nargs="+", dest="subject_list", required=True)
        group.add_argument('--workingdir', action=usage.store_directory, required=True)
        group.add_argument('--plugin', nargs="+", action=usage.store_plugin, required=True)
        group.add_argument('--name')
        
        group = parser.add_argument_group('feat.py Options')
        group.add_argument('-c', '--config', type=usage.store_filename, required=True)
        group.add_argument("--combine", action="store_true")
        group.add_argument("--fsf", action="store_true")
        group.add_argument("--feat", action="store_true")
        group.add_argument("--regress", action="store_true")
        group.add_argument("--verbose", action="store_true")
        group.add_argument("--debug", action="store_true")
        group.add_argument("--dry-run", action="store_true")
        
        return parser
    

def main(arglist):
    pp = FeatParser()
    pp.run(feat, arglist)

if __name__ == "__main__":
    main(sys.argv[1:])

