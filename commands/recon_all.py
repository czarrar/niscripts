#!/usr/bin/env python

import nipype.interfaces.io as nio # Data i/o
import nipype.interfaces.freesurfer as fs # freesurfer
import nipype.interfaces.utility as util # utility
import nipype.pipeline.engine as pe # pypeline engine
from nipype.interfaces.base import Undefined

import os, argparse, sys
import os.path as op
import numpy as np

sys.path.append(os.path.join(os.environ.get("NISCRIPTS"), "include"))
import usage
from zlogger import (LoggerError, LoggerCritical)

def recon_all(workingdir, name="recon_all_wf", **kwrds):
    
    wf = pe.Workflow(name=name)
    wf.base_dir = workingdir
    
    subinfo = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name='subinfo', 
                        iterables=('subject_id', subject_list))
    
    reconall = fs.ReconAll(ignore_exception=True)
    for k,v in kwrds:
        setattr(reconall.inputs, k, v)
    fsnode = pe.Node(interface=reconall, name="reconall")
    
    wf.connect(subinfo, 'subject_id', fsnode, 'subject_id')
    
    return wf

class ReconAllParser(usage.NiParser):
    def _create_parser(self, *args, **kwrds):
        """Create command-line interface"""
        parser = usage.NiArgumentParser(fromfile_prefix_chars="@", 
                                  argument_default=argparse.SUPPRESS,
                                  description="""Wrapper around recon-all that parallelizes processing of different subjects.""")
        parser._add_inputs = False
        parser._add_outputs = False
        
        group = parser.add_argument_group('General Options')
        group.add_argument('-s', '--subjects', nargs="+", dest="subject_list", required=True)
        group.add_argument('--workingdir', action=usage.store_directory, required=True)
        group.add_argument('--plugin', nargs="+", action=usage.store_plugin, required=True)
        group.add_argument('--name')
        
        group = parser.add_argument_group('Recon All Options')
        group.add_argument('-d', '--directive', choices=['all', 'autorecon1', 'autorecon2', 'autorecon2-cp', 'autorecon2-wm', 'autorecon2-inflate1', 'autorecon2-perhemi', 'autorecon3', 'localGI', 'qcache'], required=True)
        group.add_argument('--subjects-dir', action=usage.store_directory, required=True)
        group.add_argument('-i', '--inputs', nargs='+', dest="T1_files")
        group.add_argument('--args')
        group.add_argument('--flags')
        group.add_argument('--hemi', choices=['lh', 'rh'])
        
        return parser
    

def main(arglist):
    pp = ReconAllParser()
    pp.run(recon_all, arglist)

if __name__ == "__main__":
    main(sys.argv[1:])

