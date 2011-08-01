#!/usr/bin/env python
"""
This script uses nipype for anatomical preprocessing
"""


import numpy as np

from nipype.utils.filemanip import fname_presuffix
import nipype.interfaces.io as nio # Data i/o
import nipype.interfaces.fsl as fsl # fsl
import nipype.interfaces.afni as afni # afni
from nipype.interfaces.afni.base import Info
import nipype.interfaces.utility as util # utility
import nipype.pipeline.engine as pe # pypeline engine
import os # system functions
import argparse # command-line
import os.path as op

import sys
sys.path.append(op.join(op.dirname(op.abspath( __file__ )), "../include"))

import e_afni, misc, usage_subjects # my own extra stuff



def create_anat_preproc_workflow(subject_list, input_basedir, input_prefix, output_basedir, output_anatdir, workingdir, output_type="NIFTI_GZ", name="anatomical_preprocessing"):
    
    ######
    # Setup data source
    ######
    
    # File extension business
    ext = Info.outputtype_to_ext(output_type)
    afni.AFNICommand.set_default_outputtype(output_type)
    
    # Subject source node
    subinfo = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name='subinfo', 
                        iterables=('subject_id', subject_list))
                    
    # Say where to find input data
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['struct']), name='datasource')
    datasource.inputs.base_directory=os.path.abspath(input_basedir)
    datasource.inputs.template = os.path.join("%s", "%s%s" % (input_prefix, ext))
    datasource.inputs.template_args = dict(struct=[['subject_id']])
    
    
    ######
    # Commands
    ######
    
    # Main processes
    copy = pe.Node(interface=e_afni.Threedcopy(), name='copy')
    deoblique = pe.Node(interface=afni.ThreedWarp(deoblique=True), name='deoblique')
    reorient = pe.Node(interface=e_afni.Threedresample(orientation='RPI'), name='reorient')
    skull_strip = pe.Node(interface=e_afni.ThreedSkullStrip(), name='skull_strip')
    brain_mask = pe.Node(interface=fsl.ImageMaths(op_string='-bin'), name='brain_mask')
    
    # Pictures!
    ## of head
    slicer_head = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                            name='slicer_head')
    ## of brain mask on top of the head
    slicer_mask1 = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), name='slicer_mask1')
    slicer_mask2 = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="sagittal"), name='slicer_mask2')
    def get_overlay_args(fname):
        """Args for the overlay1 option of slicer"""
        return (fname, 1, 1)
    
    ######
    # Setup Output
    ######
    
    # Rename filenames
    headname = pe.Node(util.Rename(format_string="head", keep_ext=True), 
                        name="headname")
    brainname = pe.Node(util.Rename(format_string="brain", keep_ext=True), 
                        name="brainname")
    maskname = pe.Node(util.Rename(format_string="brain_mask", keep_ext=True), 
                        name="maskname")
    pic_headname = pe.Node(util.Rename(format_string="pic_head_axial", keep_ext=True), 
                            name="pic_headname")
    pic_maskname1 = pe.Node(util.Rename(format_string="pic_brain_mask_axial", keep_ext=True), 
                            name="pic_maskname1")
    pic_maskname2 = pe.Node(util.Rename(format_string="pic_brain_mask_sagittal", keep_ext=True), 
                            name="pic_maskname2")
    
    
    # Datasink
    ## will get: "base_directory/subject_id/output_anatdir"
    datasink = pe.Node(interface=nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = os.path.abspath(output_basedir)
    def get_substitutions(sid):
        return ('_subject_id_%s' % sid, "_subject_id_XXX")
    datasink.inputs.regexp_substitutions = ("_subject_id_XXX", output_anatdir)
    
    
    ######
    # Do the workflow!
    ######
    
    ap_pipeline = pe.Workflow(name=name)
    ap_pipeline.base_dir = os.path.abspath(workingdir)
    
    ap_pipeline.connect([(subinfo, datasource, [('subject_id', 'subject_id')]),
                        (datasource, copy, [('struct', 'in_file')]), 
                        (copy, deoblique, [('out_file', 'in_file')]), 
                        (deoblique, reorient, [('out_file', 'in_file')]),
                        (reorient, skull_strip, [('out_file', 'in_file')]),
                        (skull_strip, brain_mask, [('out_file', 'in_file')]),
                        (reorient, slicer_head, [('out_file', 'in_file')]),
                        (reorient, slicer_mask1, [('out_file', 'in_file')]),
                        (reorient, slicer_mask2, [('out_file', 'in_file')]),
                        (brain_mask, slicer_mask1, [(('out_file', get_overlay_args), 'overlay1')]),
                        (brain_mask, slicer_mask2, [(('out_file', get_overlay_args), 'overlay1')]),
                        (reorient, headname, [('out_file', 'in_file')]), 
                        (skull_strip, brainname, [('out_file', 'in_file')]), 
                        (brain_mask, maskname, [('out_file', 'in_file')]), 
                        (slicer_head, pic_headname, [('out_file', 'in_file')]),
                        (slicer_mask1, pic_maskname1, [('out_file', 'in_file')]),
                        (slicer_mask2, pic_maskname2, [('out_file', 'in_file')]),
                        (subinfo, datasink, [('subject_id', 'container'), 
                         (('subject_id', get_substitutions), 'substitutions')]), 
                        (headname, datasink, [('out_file', '@head')]),
                        (brainname, datasink, [('out_file', '@brain')]),
                        (maskname, datasink, [('out_file', '@brain_mask')]),
                        (pic_headname, datasink, [('out_file', '@pic_head')]),
                        (pic_maskname1, datasink, [('out_file', '@pic_brain_mask1')]),
                        (pic_maskname2, datasink, [('out_file', '@pic_brain_mask2')])
                        ])
    
    return ap_pipeline


def create_parser():
    """Create command-line interface"""
    parser = argparse.ArgumentParser(
        description="""
            Run preprocessing for each participant's structual images.
        """,
        parents=[usage_subjects.parent_parser]
    )
    group = parser.add_argument_group('Program Specific Options')
    group.add_argument('-i', '--input-prefix', required=True)
    group.add_argument('-o', '--output-dir', required=True, dest="output_anatdir")
    return parser


def main(arglist):
    parser = create_parser()
    args = parser.parse_args(arglist)
    kwrds = vars(args)
    plugin = kwrds.pop('plugin'); plugin_args = kwrds.pop('plugin_args')
    ap_pipeline = create_anat_preproc_workflow(**kwrds)
    ap_pipeline.run(plugin=plugin, plugin_args=plugin_args)
    return


def test_ap():
    """Testing on my computer"""
    subject_list = ["sub05676"]
    input_basedir = "/Users/zarrar/Projects/tnetworks/"
    input_prefix = "anat/mprage"
    output_basedir = "/Users/zarrar/Projects/tnetworks/output/"
    output_anatdir = "anat"
    output_workingdir = "/Users/zarrar/Projects/tnetworks/nipype/"
    ap_pipeline = create_anat_preproc_workflow(subject_list, input_basedir, input_prefix, output_basedir, output_anatdir, output_workingdir)
    ap_pipeline.run()
    ap_pipeline.write_graph()
    return

if __name__ == "__main__":
    main(sys.argv[1:])
