#!/usr/bin/env python

import nipype.interfaces.utility as util # utility
import nipype.pipeline.engine as pe # pypeline engine

class SimpleOutputConnector(object):
    """Regular Node (non Map Node version)"""
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

