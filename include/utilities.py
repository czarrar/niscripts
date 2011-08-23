#!/usr/bin/env python

import re
import nipype.interfaces.utility as util # utility
import nipype.pipeline.engine as pe # pypeline engine

__reExt = re.compile("[.]([a-zA-Z0-9]+)[.](gz|bz2)$|[.]([a-zA-Z0-9]+)$")

def sink_outputs(workflow, outputnode, datasinknode, name_prefix=""):
    if name_prefix != "" and not name_prefix.endswith("_"):
        name_prefix += "_"
    outfields = outputnode.outputs.get()
    for field in outfields:
        workflow.connect(outputnode, field, datasinknode, "@%s%s" % (name_prefix, field))

def remove_ext(fname):
    return __reExt.sub("", fname)

def get_ext(fname):
    return __reExt.search.group(fname)

class SimpleOutputConnector(object):
    """Regular Node (non Map Node version)"""
    def __init__(self, workflow, outnode):
        self.workflow = workflow
        self.outnode = outnode
        return
    
    def __call__(self, *args, **kwrds):
        return self.connect(*args, **kwrds)
    
    def connect(self, procnode, procfield, outfield, outnode=None, to_map=True, **rename_kwrds):
        default_format_string = re.sub("_pic$", "", outfield)
        rename_kwrds.setdefault('format_string', default_format_string)
        rename_kwrds.setdefault('keep_ext', True)
        rename = util.Rename(**rename_kwrds)
        if outnode is None:
            outnode = self.outnode
        renamenode = pe.Node(interface=rename, 
                             name="rename_%s" % outfield)            
        self.workflow.connect([
            (procnode, renamenode, [(procfield, "in_file")]), 
            (renamenode, outnode, [("out_file", outfield)])
        ])
        return renamenode
    



class OutputConnector(object):
    def __init__(self, workflow, outnode):
        self.workflow = workflow
        self.outnode = outnode
        return
    
    def __call__(self, *args, **kwrds):
        return self.connect(*args, **kwrds)
    
    def connect(self, procnode, procfield, outfield, outnode=None, to_map=True, **rename_kwrds):
        rename_kwrds.setdefault('format_string', outfield)
        rename_kwrds.setdefault('keep_ext', True)
        rename = util.Rename(**rename_kwrds)
        if outnode is None:
            outnode = self.outnode
        if to_map:
            renamenode = pe.MapNode(interface=rename, iterfield=["in_file"], 
                                    name="rename_%s" % outfield)
        else:
            renamenode = pe.Node(interface=rename, 
                                 name="rename_%s" % outfield)            
        self.workflow.connect([
            (procnode, renamenode, [(procfield, "in_file")]), 
            (renamenode, outnode, [("out_file", outfield)])
        ])
        return renamenode
    

# borrowed from Fluid_NiPype
def get_mapnode_substitutions(workflow, output_node, nruns, unmap=[]):
    import networkx as nx
    from nipype.pipeline.engine import MapNode
    substitutions = []
    unmapfields = [ "rename_%s" % outfield for outfield in unmap ]
    mapnodes = [
        e[0].name for e in nx.edges(workflow._graph) \
            if e[1] is output_node and isinstance(e[0],MapNode) and e[0].name not in unmapfields
    ]
    unmapnodes = [
        e[0].name for e in nx.edges(workflow._graph) \
            if e[1] is output_node and isinstance(e[0],MapNode) and e[0].name in unmapfields
    ]
    for r in range(nruns):
        for node in mapnodes:
            substitutions.append(("_%s%d"%(node, r), "run_%d"%(r+1)))
        for node in unmapnodes:
            substitutions.append(("_%s%d"%(node, r), ""))
    return substitutions

def run_freesurfer_fun(subjects_dir, subject_id, directive, t1_files=[]):
    import os
    import os.path as op
    import nipype.interfaces.freesurfer as fs
    from execute import Process
    from glob import glob
    
    subdir = op.join(subjects_dir, subject_id)
    recon = fs.ReconAll(directive='autorecon1', subjects_dir=subjects_dir, 
                        subject_id=subject_id)
    
    rawfile = op.join(subdir, "mri", "orig", "001.mgz")
    if t1_files:
        if op.isfile(rawfile):
            print 'T1 files specified but not to be included since they already exist in output'
        else:
            recon.inputs.t1_files = t1_files
    
    to_run = False
    filechoices = {
        'autorecon1': op.join(subdir, "mri", "brainmask.mgz"),
        'autorecon2': op.join(subdir, "surf", "?h.inflated"),
        'autorecon3': op.join(subdir, "mri", "aparc+aseg.mgz")
    }
    
    try:
        fpath = filechoices[directive]
    except KeyError:
        raise Exception("Unrecognized directive %s" % directive)
    
    ran = False
    gpath = glob(fpath)
    if len(gpath) > 0:
        print 'Freesurfer output for directive %s already exists' % directive
    else:
        recon.inputs.directive = directive
        p = Process(recon.cmdline, to_print=True)
        if p.retcode != 0:
            raise Exception("Error running '%s': \n%s" % (recon.cmdline, p.stderr))
        ran = True
    
    return (subject_id, ran)


run_freesurfer = util.Function(input_names=["subjects_dir", "subject_id", "directive", 
                                            "t1_files"], 
                               output_names=["subject_id", "ran"], 
                               function=run_freesurfer_fun)
