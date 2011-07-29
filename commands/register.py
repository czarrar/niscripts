#!/usr/bin/env python
"""
This script uses nipype for normalization
"""

import sys
sys.path.append('/Users/zarrar/Code/nipype') # replace with dynamic path
sys.path.append('/Users/zarrar/Dropbox/code/niscripts/include') # replace with dynamic path

import numpy as np

from nipype.interfaces.base import Bunch
from nipype.utils.filemanip import fname_presuffix
import nipype.interfaces.io as nio # Data i/o
import nipype.interfaces.fsl as fsl # fsl
import nipype.interfaces.afni as afni # afni
from nipype.interfaces.afni.base import Info
import nipype.interfaces.utility as util # utility
import nipype.pipeline.engine as pe # pypeline engine
import os # system functions
import argparse # command-line
import re

import e_afni, misc, usage_subjects # my own extra stuff
from utilities import *

tolist = lambda x: [x]

def reg_pics(in_file, ref_file):
    COMMANDS[6]="slicer fnirt_highres2standard standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png"
    COMMANDS[7]="pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png fnirt_highres2standard1.png"
    COMMANDS[8]="slicer standard fnirt_highres2standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png"
    COMMANDS[9]="pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png fnirt_highres2standard2.png"
    COMMANDS[10]="pngappend fnirt_highres2standard1.png - fnirt_highres2standard2.png fnirt_highres2standard.png"
    

class RegOutputConnector(object):
    def __init__(self, workflow, outnode, innode):
        self.innode = innode
        self.workflow = workflow
        self.outnode = outnode
        return
    
    def __call__(self, *args, **kwrds):
        return self.connect(*args, **kwrds)
    
    def connect(self, procnode, procfield, outfield, **rename_kwrds):
        ftype = re.sub("_.*", "", outfield)
        if ftype == "in2ref" or ftype == "ref2in":
            rename_kwrds.setdefault('format_string', "%(a)s2%(b)s")
        else:
            rename_kwrds.setdefault('format_string', outfield)
        
        rename_kwrds.setdefault('keep_ext', True)
        rename = util.Rename(**rename_kwrds)
        rnode = pe.Node(interface=rename, iterfield=["in_file"], 
                                name="rename_%s" % outfield)
        self.workflow.connect([
            (procnode, rnode, [(procfield, "in_file")]), 
            (rnode, self.outnode, [("out_file", outfield)])
        ])
        
        if ftype == "in2ref":
            self.workflow.connect(self.innode, 'in_prefix', rnode, 'a')
            self.workflow.connect(self.innode, 'ref_prefix', rnode, 'b')
        elif ftype == "ref2in"
            self.workflow.connect(self.innode, 'ref_prefix', rnode, 'a')
            self.workflow.connect(self.innode, 'in_prefix', rnode, 'b')
        
        return rnode


# applywarp
('nn' or 'trilinear' or 'sinc' or 'spline')

def create_nonlin_reg_workflow(
    name = "nonlinear_registration", 
    regtype = 'highres2standard'):
    
    # check reg_type
    reg_opts = ['highres2standard', 'manual']
    if reg_type not in reg_opts:
        raise Exception("reg_type must be one of: %s" % ", ".join(reg_opts))
    
    
    #####
    # Setup workflow
    #####
    
    normalize = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        # required
        "in_file",      # head
        "ref_file",     # standard head
        "affine_file",  # mat: highres -> standard
        "refmask_file", # standard brain mask
        "config_file",  # file specifying command-line parameters
        "in_prefix",
        "out_prefix",
        # internal
        "wtype" # lin, nonlin, lin_linker, or nonlin_linker
    ]
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields), 
                            name="inputspec")
    
    # set workflow type
    inputnode.inputs.wtype = 'nonlin'
    
    # regtype
    if reg_type != 'manual':
        re_io = re.search("(?P<in>\w+)2(?P<out>\w+)", reg_type).groupdict()
        inputnode.inputs.in_prefix = re_io['in']
        inputnode.inputs.ref_prefix = re_io['out']
    
    #####
    # Setup output node
    #####
    
    # Outputs
    output_fields = [
        "log_file",
        "in2ref_warped",
        "fieldcoeff_file",
        "jacobian_file",
        #"in2ref_pic"
    ]
    outputnode = pe.Node(util.IdentityInterface(fields=output_fields),
                        name="outputspec")
    
    # Rename output filenames
    renamer = RegOutputConnector(normalize, outputnode, inputnode)
    
    
    #####
    # Commands
    #####
    
    fnirt = pe.Node(fsl.FNIRT(fieldcoeff_file=fieldcoeff, jacobian_file=jacobian),
                    name="fnirt")
    normalize.connect([
        (inputnode, fnirt, [('in_file', 'in_file'), 
                            ('ref_file', 'ref'), 
                            ('affine_file', 'affine_file'), 
                            ('refmask_file', 'refmask_file'), 
                            ('config_file', 'config_file')])
    ])
    renamer(fnirt, 'log_file', 'log_file', format_string="fnirt_log")
    renamer(fnirt, 'warped_file', 'in2ref_warped')
    renamer(fnirt, 'fieldcoef_file', 'fieldcoeff_file', format_string"%(a)s2%(b)s_warp")
    renamer(fnirt, 'jacobian_file', 'jacobian_file', format_string"%(a)s2%(b)s_jac")
        
    # pics (using slicer)
    ## input over ref
    reg_slicer_a = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                            name='reg_slicer_a')
    normalize.connect([
        (inputnode, reg_slicer_a, [('ref_file', 'in_file')]),
        (fnirt, reg_slicer_a, [('out_file', 'edge_overlay')]),

    ])
    ## ref over input
    reg_slicer_b = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                            name='reg_slicer_b')
    normalize.connect([
        (fnirt, reg_slicer_b, [('warped_file', 'in_file')]),
        (inputnode, reg_slicer_b, [('ref_file', 'edge_overlay')])
    ])
    # combine
    #concat_pics = pe.Node(interface=util.Merge(2, axis="hstack"), 
    #                name="concat_pics")
    #normalize.connect([
    #    (reg_slicer_a, concat_pics, [(('out_file', tolist), 'in1')]),
    #    (reg_slicer_b, concat_pics, [(('out_file', tolist), 'in2')])
    #])
    #reg_pngappend = pe.Node(interface=fsl.Pngappend(op_string="%s - %s"), 
    #                        name="reg_pngappend")
    #normalize.connect(concat_pics, 'out', reg_pngappend, 'in_files')
    #renamer(reg_pngappend, 'out_file', 'in2ref_png')
    
    return normalize


# have search 
# search_type: 'nada', 'normal', 'full', 'manual'
# reg_type: func2init, func2highres, init2highres, highres2standard, None
def create_lin_reg_workflow(
    name="linear_registration", 
    reg_type = 'manual', 
    search_type = 'normal', 
    init_matrix = False):
    """Generate workflow for linear registration"""
    
    # check reg_type
    reg_opts = ['func2coplanar', 'func2highres', 'coplanar2highres', 'highres2standard']
    if reg_type not in reg_opts:
        raise Exception("reg_type must be one of: %s" % ", ".join(reg_opts))
    
    # check search_type
    search_opts = ['nada', 'normal', 'full', 'manual']
    if search_type not in search_opts:
        raise Exception("search_type must be one of: %s" % ", ".join(search_opts))
    
    
    #####
    # Setup workflow
    #####
    
    normalize = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        # required
        "in_file",
        "ref_file",
        # required if reg_type not specified
        "in_prefix",    # e.g. highres
        "ref_prefix",   # e.g. standard
        "dof",
        # optional
        "search_x",
        "search_y",
        "search_z",
        "cost",
        "interp",
        # internal
        "wtype" # lin, nonlin, lin_linker, or nonlin_linker
    ]
    if init_matrix:
        input_fields.append("init_matrix")
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields), 
                            name="inputspec")
    
    # set defaults
    inputnode.inputs.cost = 'corratio'
    inputnode.inputs.interp = 'trilinear'
    if search_type != "manual":
        sw_opts = {
            'nada':     [0, 0], 
            'normal':   [-90,90], 
            'full':     [-180,180],   
        }
        inputnode.inputs.search_x = sw_opts[search_type]
        inputnode.inputs.search_y = sw_opts[search_type]
        inputnode.inputs.search_z = sw_opts[search_type]
    
    # set other defaults based on reg_type
    if reg_type != 'manual':
        dof_opts = dict(zip(reg_opts, [6,6,6,12]))
        inputnode.inputs.dof = dof_opts[reg_type]
        re_io = re.search("(?P<in>\w+)2(?P<out>\w+)", reg_type).groupdict()
        inputnode.inputs.in_prefix = re_io['in']
        inputnode.inputs.ref_prefix = re_io['out']
        if re_io['in'] == "func":
            inputnode.inputs.cost = 'normmi'    # although do a corratio first!
    
    # set workflow type
    inputnode.inputs.wtype = 'lin'
    
    
    #####
    # Setup output node
    #####
    
    # Outputs
    output_fields = [
        "in2ref",
        #"ref2in",
        "in2ref_mat",
        "ref2in_mat",
        #"in2ref_png",
    ]
    outputnode = pe.Node(util.IdentityInterface(fields=output_fields),
                        name="outputspec")
    
    # Rename output filenames
    renamer = RegOutputConnector(normalize, outputnode, inputnode)
    #renamer = OutputConnector(normalize, outputnode)
    
    #def easier_renamer(procnode, procfield, outfield, **kwrds):
    #    rnode = renamer.connect(procnode, procfield, outfield, 
    #                            format_string="%(a)s2%(b)s", **kwrds)
    #    ftype = re.sub("_.*", outfield)
    #    if ftype == "in2ref":
    #        normalize.connect(outputnode, 'in_prefix', rnode, 'a')
    #        normalize.connect(outputnode, 'ref_prefix', rnode, 'b')
    #    else:
    #        normalize.connect(outputnode, 'ref_prefix', rnode, 'a')
    #        normalize.connect(outputnode, 'in_prefix', rnode, 'b')
    #    
    #    return rnode
    
    #####
    # Commands
    #####
    
    if reg_type and re_io['in'] == 'func':  # do a 2-stage registration
        flirt_init = pe.Node(fsl.FLIRT(cost='corratio', cost_func='corratio'), 
                            name="flirt_init")
        normalize.connect([
            (inputnode, flirt_init, [('in_file', 'in_file')]),
            (inputnode, flirt_init, [('ref_file', 'reference')]),
            (inputnode, flirt_init, [('dof', 'dof')]),
            (inputnode, flirt_init, [('search_x', 'searchr_x')]),
            (inputnode, flirt_init, [('search_y', 'searchr_y')]),
            (inputnode, flirt_init, [('search_z', 'searchr_z')])
        ])
        if init_matrix:
            normalize.connect(inputnode, 'init_matrix', flirt_init, 'in_matrix_file')
    
    flirt = pe.Node(fsl.FLIRT(), name="flirt")
    normalize.connect([
        (inputnode, flirt, [('in_file', 'in_file')]),
        (inputnode, flirt, [('ref_file', 'reference')]),
        (inputnode, flirt, [('cost', 'cost')]),
        (inputnode, flirt, [('cost', 'cost_func')]),
        (inputnode, flirt, [('dof', 'dof')]),
        (inputnode, flirt, [('search_x', 'searchr_x')]),
        (inputnode, flirt, [('search_y', 'searchr_y')]),
        (inputnode, flirt, [('search_z', 'searchr_z')]),
        (inputnode, flirt, [('interp', 'interp')])
    ])
    renamer(flirt, 'out_file', 'in2ref')
    renamer(flirt, 'out_matrix_file', 'in2ref_mat')
        
    if reg_type and re_io['in'] == 'func':
        normalize.connect(flirt_init, 'out_matrix_file', flirt, 'in_matrix_file')
    elif init_matrix:
        normalize.connect(inputnode, 'init_matrix', flirt, 'in_matrix_file')
    
    # convert_xfm (to get ref2in_mat)
    invert_mat = pe.Node(fsl.ConvertXFM(invert_xfm=True), name="invert_mat")
    normalize.connect(flirt, 'out_matrix_file', invert_mat, 'in_file')
    renamer(flirt, 'out_matrix_file', 'ref2in_mat')
        
    ## applyxfm (to get ref2in)
    #apply_invert_mat = pe.Node(fsl.ApplyXfm(apply_xfm=True), name="apply_invert_mat")
    #normalize.connect([
    #    (inputnode, apply_invert_mat, [('ref_file', 'in_file')]),
    #    (inputnode, apply_invert_mat, [('in_file', 'reference')]),
    #    (invert_mat, apply_invert_mat, [('out_file', 'in_matrix_file')]),
    #    (inputnode, apply_invert_mat, [('interp', 'interp')])
    #])
    #renamer(apply_invert_mat, 'out_matrix_file', 'ref2in')
    
    # pics (using slicer)
    ## input over ref
    reg_slicer_a = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                            name='reg_slicer_a')
    normalize.connect([
        (inputnode, reg_slicer_a, [('ref_file', 'in_file')]),
        (flirt, reg_slicer_a, [('out_file', 'edge_overlay')]),
        
    ])
    ## ref over input
    reg_slicer_b = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                            name='reg_slicer_b')
    normalize.connect([
        (flirt, reg_slicer_b, [('out_file', 'in_file')]),
        (inputnode, reg_slicer_b, [('ref_file', 'edge_overlay')])
    ])
    # combine
    #concat_pics = pe.Node(interface=util.Merge(2, axis="hstack"), 
    #                name="concat_pics")
    #normalize.connect([
    #    (reg_slicer_a, concat_pics, [(('out_file', tolist), 'in1')]),
    #    (reg_slicer_b, concat_pics, [(('out_file', tolist), 'in2')])
    #])
    #reg_pngappend = pe.Node(interface=fsl.Pngappend(op_string="%s - %s"), 
    #                        name="reg_pngappend")
    #normalize.connect(concat_pics, 'out', reg_pngappend, 'in_files')
    #renamer(reg_pngappend, 'out_file', 'in2ref_png')
    
    return normalize


def create_link_lin_reg_workflow(
    name="link_linear_registrations",
    in2x_reg_type="manual", 
    x2ref_reg_type="manual"):
    """does stuff"""
    
    # check reg_type
    reg_opts = ['func2coplanar', 'func2highres', 'coplanar2highres', 'highres2standard']
    if in2x_reg_type not in reg_opts or x2ref_reg_type not in reg_opts:
        raise Exception("in2x_reg_type or x2ref_reg_type must be one of: %s" % ", ".join(reg_opts))
    
    
    #####
    # Setup workflow
    #####
    
    linker = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        # required
        "in_file",
        "in2x_mat_file",
        "x2ref_mat_file",
        "ref_file",
        # required if reg_types not specified
        "in_prefix",    # e.g. highres
        "ref_prefix",    # e.g. standard
        # optional
        "interp",
        # internal
        "wtype" # lin, nonlin, lin_linker, or nonlin_linker
    ]
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields), 
                            name="inputspec")
    
    # set defaults
    inputnode.inputs.interp = 'trilinear'
    inputnode.inputs.wtype = 'lin_linker'
    
    # set in and out prefix
    if in2x_reg_type != 'manual':
        re_io = re.search("(?P<in>\w+)2(?P<x>\w+)", in2x_reg_type).groupdict()
        inputnode.inputs.in_prefix = re_io['in']
    if x2ref_reg_type != 'manual':
        re_io = re.search("(?P<x>\w+)2(?P<ref>\w+)", x2ref_reg_type).groupdict()
        inputnode.inputs.ref_prefix = re_io['ref']
    
    
    #####
    # Setup output node
    #####
    
    # Outputs
    output_fields = [
        "in2ref",
        #"ref2in",
        "in2ref_mat",
        "ref2in_mat",
        #"in2ref_png",
    ]
    outputnode = pe.Node(util.IdentityInterface(fields=output_fields),
                        name="outputspec")
    
    # Rename output filenames
    renamer = RegOutputConnector(linker, outputnode, inputnode)
    
    # convert_xfm
    concat_mat = pe.Node(fsl.ConvertXFM(concat_xfm=True), name="concat_mat")
    linker.connect(inputnode, 'in2x_mat_file', concat_mat, 'in_file')
    linker.connect(inputnode, 'x2ref_mat_file', concat_mat, 'in_file2')
    renamer(concat_mat, 'out_file', 'in2ref_mat')
        
    # convert_xfm
    invert_mat = pe.Node(fsl.ConvertXFM(invert_xfm=True), name="invert_mat")
    linker.connect(concat_mat, 'out_file', invert_mat, 'in_file')
    renamer(invert_mat, 'out_file', 'ref2in_mat')
    
    # applyxfm
    apply_mat = pe.Node(fsl.ApplyXfm(apply_xfm=True), name="apply_mat")
    linker.connect([
        (inputnode, apply_mat, [('in_file', 'in_file')]),
        (inputnode, apply_mat, [('ref_file', 'reference')]),
        (concat_mat, apply_mat, [('out_file', 'in_matrix_file')]),
        (inputnode, apply_mat, [('interp', 'interp')])
    ])
    renamer(apply_mat, 'out_file', 'in2ref')
    
    # pics (using slicer)
    ## input over ref
    reg_slicer_a = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                            name='reg_slicer_a')
    linker.connect([
        (inputnode, reg_slicer_a, [('ref_file', 'in_file')]),
        (apply_mat, reg_slicer_a, [('out_file', 'edge_overlay')]),
        
    ])
    ## ref over input
    reg_slicer_b = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                            name='reg_slicer_b')
    linker.connect([
        (apply_mat, reg_slicer_b, [('out_file', 'in_file')]),
        (inputnode, reg_slicer_b, [('ref_file', 'edge_overlay')])
    ])
    ### combine
    #concat_pics = pe.Node(interface=util.Merge(2, axis="hstack"), 
    #                name="concat_pics")
    #linker.connect([
    #    (reg_slicer_a, concat_pics, [(('out_file', tolist), 'in1')]),
    #    (reg_slicer_b, concat_pics, [(('out_file', tolist), 'in2')])
    #])
    #reg_pngappend = pe.Node(interface=fsl.Pngappend(op_string="%s - %s"), 
    #                        name="reg_pngappend")
    #linker.connect(concat_pics, 'out', reg_pngappend, 'in_files')
    #renamer(reg_pngappend, 'out_file', 'in2ref_png')
    
    return linker


def create_func2standard_workflow(
    name="func2standard",
    coplanar = False, 
    search_type = 'normal'):
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        "func",
        "coplanar", # optional
        "highres",
        "standard",
        "interp"    # optional
    ]
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields, 
                                                         mandatory_inputs=False), 
                            name="inputspec")
    
    # set defaults
    inputnode.inputs.interp = 'trilinear'
    
    
    #####
    # Setup workflow
    #####
    
    normalize = pe.Workflow(name=name)
    workflows = Bunch(func=[], highres=[])
    
    
    #####
    # Commands
    #####
    
    # highres2standard
    name = "highres2standard"
    highres2standard = create_lin_reg_workflow(name = name, 
                                               reg_type = name, 
                                               search_type = search_type)
    workflows.highres.append(highres2standard)
    normalize.connect([
        (inputnode, highres2standard, [('highres', 'inputspec.in_file'), 
                                       ('standard', 'inputspec.ref_file')])
    ])
    
    if coplanar:
        # func2coplanar
        name1 = "func2coplanar"
        func2coplanar = create_lin_reg_workflow(name = name1, 
                                                reg_type = name1, 
                                                search_type = search_type)
        workflows.func.append(func2coplanar)
        # coplanar2highres
        name2 = "coplanar2highres"
        coplanar2highres = create_lin_reg_workflow(name = name2,
                                                   reg_type = name2, 
                                                   search_type = search_type)
        workflows.func.append(coplanar2highres)
        # func2highres
        name3 = "func2highres"
        func2highres = create_link_lin_reg_workflow(name = name3, 
                                                    in2x_reg_type = name1, 
                                                    x2ref_reg_type = name2)
        workflows.func.append(func2highres)
        # link them
        normalize.connect([
            (inputnode, func2coplanar, [('func', 'inputspec.in_file'), 
                                        ('coplanar', 'inputspec.ref_file')]),
            (inputnode, coplanar2highres, [('coplanar', 'inputspec.in_file'), 
                                           ('highres', 'inputspec.ref_file')]),
            (inputnode, func2highres, [('func', 'inputspec.in_file'), 
                                       ('highres', 'inputspec.ref_file')]),
            (func2coplanar, func2highres, [('outputspec.in2ref_mat', 
                                            'inputspec.in2x_mat_file')]), 
            (coplanar2highres, func2highres, [('outputspec.in2ref_mat', 
                                               'inputspec.x2ref_mat_file')])
        ])
    else:
        # func2highres 
        name = "func2highres"
        func2highres = create_lin_reg_workflow(name = name, 
                                               reg_type = name, 
                                               search_type = search_type)
        workflows.func.append(func2highres)
        normalize.connect([
            (inputnode, func2highres, [('func', 'inputspec.in_file'), 
                                       ('highres', 'inputspec.ref_file')])
        ])
    
    # func2standard
    name = "func2standard"
    func2standard = create_link_lin_reg_workflow(name = name, 
                                                 in2x_reg_type = "func2highres", 
                                                 x2ref_reg_type = "highres2standard")
    workflows.func.append(func2standard)
    normalize.connect([
        (inputnode, func2standard, [('func', 'inputspec.in_file'), 
                                    ('standard', 'inputspec.ref_file')]),
        (func2highres, func2standard, [('outputspec.in2ref_mat', 'inputspec.in2x_mat_file')]),
        (highres2standard, func2standard, [('outputspec.in2ref_mat', 'inputspec.x2ref_mat_file')])
    ])
    
    # link interp input
    for wf in workflows.func:
        normalize.connect(inputnode, 'interp', wf, 'inputspec.interp')    
    for wf in workflows.highres:
        normalize.connect(inputnode, 'interp', wf, 'inputspec.interp')    
    
    return (normalize, workflows)


# inputs (bunch): func, coplanar, highres
# outputs (bunch): func, coplanar, highres
def register( 
    subject_list, 
    input_basedir, inputs, standard, 
    output_basedir, outputs, workingdir, 
    interp, search, output_type, 
    name="linear_registration_func2standard"):
    
    if isinstance(inputs, dict):
        inputs = Bunch(**inputs)
    if isinstance(outputs, dict):
        outputs = Bunch(**outputs)
    # todo: check elements in inputs and outputs
    
    #####
    # Setup workflow
    #####
    
    if hasattr(inputs, "coplanar"):
        have_coplanar = True
    else:
        have_coplanar = False
    regproc, workflows = create_func2standard_workflow(name, have_coplanar, search)
    regproc.base_dir = workingdir
    
    # get input / set certain inputs
    inputnode = regproc.get_node("inputspec")
    inputnode.inputs.standard = standard
    inputnode.inputs.interp = interp
    
    
    ######
    # Setup data source
    ######
    
    # File extension
    ext = Info.outputtype_to_ext(output_type)
    afni.AFNICommand.set_default_outputtype(output_type)
    
    # Subject ID Node
    subinfo = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name='subinfo', 
                        iterables=('subject_id', subject_list))
    
    # Location of input data
    outfields = ['func', 'highres']
    if have_coplanar:
        outfields.append('coplanar')
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], 
                                                    outfields=outfields), 
                         name='datasource')
    datasource.inputs.base_directory=os.path.abspath(input_basedir)
    datasource.inputs.template = "*"
    datasource.inputs.field_template = dict(
        func = os.path.join("%s", inputs.func),
        highres = os.path.join("%s", inputs.highres)
    )
    datasource.inputs.template_args = dict(
                                        func = [['subject_id']],
                                        highres = [['subject_id']]
                                      )
    if have_coplanar:
        datasource.inputs.field_template['coplanar'] = os.path.join("%s", "%s%s" % (inputs.coplanar, ext))
        datasource.inputs.template_args['coplanar'] = [['subject_id']]
    
    
    # Link inputs
    regproc.connect([
        (subinfo, datasource, [('subject_id', 'subject_id')]),
        (datasource, inputnode, [('func', 'func'), ('highres', 'highres')])
    ])
    if have_coplanar:
        regproc.connect(datasource, 'coplanar', inputnode, 'coplanar')
    
    
    ######
    # Setup data sink
    ######
    
    # Datasinks
    datasinks = dict.fromkeys(vars(workflows).keys())
    ## will get: "base_directory/subject_id/output_dir"
    for k in datasinks:
        datasinks[k] = pe.Node(interface=nio.DataSink(), name='datasink_%s' % k)
        datasinks[k].inputs.base_directory = os.path.abspath(output_basedir)
        # set container to subject_id
        regproc.connect(subinfo, 'subject_id', datasinks[k], 'container')
        # replace subject_id stuff with output directory and reg
        datasinks[k].inputs.regexp_substitutions = (r"_subject_id_(\w|\d)+", getattr(outputs, k))
    
    # connect output to datasinks
    for k, wfs in vars(workflows).iteritems():
        for wf in wfs:
            outnode = wf.get_node("outputspec")
            outfields = outnode.outputs.get()
            for outfield in outfields:
                regproc.connect(outnode, outfield, datasinks[k], 
                                "@%s_%s" % (wf.name, outfield))
    
    return regproc


class store_func_highres(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        k = option_string.replace("-", "")
        if not hasattr(namespace, 'inputs'):
            namespace.inputs = {}
        if not hasattr(namespace, 'outputs'):
            namespace.outputs = {}
        namespace.inputs[k] = values[0]
        namespace.outputs[k] = values[1]
        return
    

class store_coplanar(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        k = option_string.replace("-", "")
        if not hasattr(namespace, 'inputs'):
            namespace.inputs = {}
        namespace.inputs[k] = value
    


#class store_regdir(argparse.Action):
#    def __call__(self, parser, namespace, value, option_string=None):
#        outputs = getattr(namespace, 'outputs', {})
#        outputs[option_string] = value
#        setattr(namespace, 'outputs', outputs)
#    
#

def create_parser():
    parser = argparse.ArgumentParser(
        description="""
            Register/normalize functional image to standard space
        """,
        parents=[usage_subjects.parent_parser]
    )
    group = parser.add_argument_group("Registration Options")
    group.add_argument("--interp", choices=["trilinear", "nearestneighbour", "sinc"], default="trilinear")
    group.add_argument("--search", choices=["nada", "normal", "full"], default="normal")
    group.add_argument("--func", nargs=2, action=store_func_highres, default=argparse.SUPPRESS, required=True)
    group.add_argument("--coplanar", action=store_coplanar, default=argparse.SUPPRESS)
    group.add_argument("--highres", nargs=2, action=store_func_highres, default=argparse.SUPPRESS, required=True)
    group.add_argument("--standard", default=fsl.Info.standard_image("MNI152_T1_2mm_brain.nii.gz"))
    
    return parser

def main(arglist):
    parser = create_parser()
    args = parser.parse_args(arglist)
    kwrds = vars(args)
    reg_pipe = register(**kwrds)
    reg_pipe.run()
    reg_pipe.write_graph()
    return

def test_p():
    arglist = "-s tb3417 -b /Users/zarrar/Projects/tnetworks/output --workingdir /Users/zarrar/Projects/tnetworks/tmp --func func/func_ref.nii.gz func/reg --highres highres/brain.nii.gz highres/reg".split()
    parser = create_parser()
    args = parser.parse_args(arglist)
    kwrds = vars(args)
    reg_pipe = register(**kwrds)
    return reg_pipe

if __name__ == "__main__":
    main(sys.argv[1:])

