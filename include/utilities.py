#!/usr/bin/env python

import nipype.interfaces.utility as util # utility
import nipype.pipeline.engine as pe # pypeline engine

class OutputConnector(object):
    def __init__(self, workflow, outnode):
        self.workflow = workflow
        self.outnode = outnode
        return
    
    def __call__(self, *args, **kwrds):
        return self.connect(*args, **kwrds)
    
    def connect(self, procnode, procfield, outfield, outnode=None, **rename_kwrds):
        rename_kwrds.setdefault('format_string', outfield)
        rename_kwrds.setdefault('keep_ext', True)
        rename = util.Rename(**rename_kwrds)
        if outnode is None:
            outnode = self.outnode
        renamenode = pe.MapNode(interface=rename, iterfield=["in_file"], 
                                name="rename_%s" % outfield)
        self.workflow.connect([
            (procnode, renamenode, [(procfield, "in_file")]), 
            (renamenode, outnode, [("out_file", outfield)])
        ])
        return renamenode
    

# borrowed from Fluid_NiPype
def get_mapnode_substitutions(workflow, output_node, runs, unmap=[]):
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
    for r in runs:
        for node in mapnodes:
            substitutions.append(("_%s%d"%(node, r-1), "run_%d"%(r)))
        for node in unmapnodes:
            substitutions.append(("_%s%d"%(node, r-1), ""))
    return substitutions

