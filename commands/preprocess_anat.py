#!/usr/bin/env python
"""
This script uses nipype for anatomical preprocessing
"""

import numpy as np

from nipype.utils.filemanip import fname_presuffix
import nipype.interfaces.io as nio # Data i/o
import nipype.interfaces.fsl as fsl # fsl
import nipype.interfaces.afni as afni # afni
import nipype.interfaces.freesurfer as fs # afni
from nipype.interfaces.afni.base import Info
import nipype.interfaces.utility as util # utility
import nipype.pipeline.engine as pe # pypeline engine
import os # system functions
import argparse # command-line
import os.path as op

import sys
sys.path.append(op.join(os.environ.get("NISCRIPTS"), "include"))

import e_afni, misc, usage # my own extra stuff
from utilities import SimpleOutputConnector, run_freesurfer
from execute import Process


def get_overlay_args(fname):
    """Args for the overlay1 option of slicer"""
    return (fname, 1, 1)

def anatomical_preprocessing(
    subject_list, 
    inputs, outputs, workingdir, output_type,
    run, orientation, freesurfer_dir=None, name="anatomical_preprocessing"):
    
    #####
    # Setup workflow
    #####
    
    if run == "bet":
        preproc = create_ap_bet_workflow(name=name + "_withbet")
        # get input / set certain inputs
        inputnode = preproc.get_node("inputspec")
    elif run == "freesurfer":
        preproc = create_ap_freesurfer_workflow(freesurfer_dir, name=name + "_withfreesurfer")
        # get input / set certain inputs
        inputnode = preproc.get_node("inputspec")
        inputnode.inputs.orientation = orientation
        inputnode.inputs.freesurfer_dir = freesurfer_dir
    preproc.base_dir = workingdir
    
    ######
    # Setup data source
    ######
    
    # Subject source node
    subinfo = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name='subinfo', 
                        iterables=('subject_id', subject_list))
    
    # Say where to find input data
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['struct']), name='datasource')
    datasource.inputs.base_directory=os.path.abspath(inputs.basedir)
    datasource.inputs.template = os.path.join("%s", inputs.struct)
    datasource.inputs.template_args = dict(struct=[['subject_id']])    
    
    # Link inputs
    preproc.connect([
        (subinfo, datasource, [('subject_id', 'subject_id')]),
        (datasource, inputnode, [('struct', 'struct')])
    ])
    if run == "freesurfer":
        preproc.connect([
            (subinfo, inputnode, [('subject_id', 'subject_id')])
        ])
    
    
    ###### 
    # Setup data sink
    ######
    
    # Datasink
    datasink = pe.Node(interface=nio.DataSink(), name='datasink')
    
    # set base directory
    datasink.inputs.base_directory = os.path.abspath(outputs.basedir)
    
    # replace subject_id stuff with anat
    datasink.inputs.regexp_substitutions = (r"_subject_id_(\w|\d)+", outputs.struct)
    
    # subject_id
    preproc.connect(subinfo, 'subject_id', datasink, 'container')
    
    ## connect
    outputnode = preproc.get_node("outputspec")
    output_fields = outputnode.outputs.get()
    for field in output_fields:
        preproc.connect(outputnode, field, datasink, "@%s" % field)
    
    return preproc


def create_ap_bet_workflow(name="preproc_anat_bet"):
    
    #####
    # Setup workflow
    #####
    
    preproc = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        "struct",
        "orientation"
    ]
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields), 
                            name="inputspec")
    
    # defaults
    inputnode.inputs.orientation = "RPI"    # fsl's standard brain
    
    
    #####
    # Setup output node
    #####
    
    # Outputs
    output_fields = [
        "orig",
        "head",
        "brain",
        "brain_mask",
        "head_axial_pic",
        "brain_mask_axial_pic",
        "brain_mask_sagittal_pic"        
    ]
    outputnode = pe.Node(util.IdentityInterface(fields=output_fields),
                        name="outputspec")
    
    # renaming
    renamer = SimpleOutputConnector(preproc, outputnode)
    
    # Deoblique (so image plays well with AFNI)
    deoblique = pe.Node(interface=afni.ThreedWarp(deoblique=True), name='deoblique')
    preproc.connect(inputnode, 'struct', deoblique, 'in_file')
    renamer(inputnode, 'struct', 'orig')
    
    # Reorient (so image plays well with registration)
    reorient = pe.Node(interface=afni.Threedresample(), name='reorient')
    preproc.connect(deoblique, 'out_file', reorient, 'in_file')
    preproc.connect(inputnode, 'orientation', reorient, 'orientation')
    renamer(reorient, 'out_file', 'head')
    
    # Skull Strip
    skull_strip = pe.Node(interface=fsl.BET(mask=True), name='skull_strip')
    preproc.connect(reorient, 'out_file', skull_strip, 'in_file')
    renamer(skull_strip, 'out_file', 'brain')
    renamer(skull_strip, 'mask_file', 'brain_mask')
    
    # Pic of head
    slicer_head = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                            name='slicer_head')
    preproc.connect(reorient, 'out_file', slicer_head, 'in_file')
    renamer(slicer_head, 'out_file', 'head_axial_pic')
    
    # Pic of brain mask overlaid on head (axial)
    slicer_mask1 = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), name='slicer_mask1')
    preproc.connect(reorient, 'out_file', slicer_mask1, 'in_file')
    preproc.connect(skull_strip, ('mask_file', get_overlay_args), slicer_mask1, 'overlay1')
    renamer(slicer_mask1, 'out_file', 'brain_mask_axial_pic')
    
    # Pic of brain mask overlaid on head (sagittal)
    slicer_mask2 = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="sagittal"), name='slicer_mask2')
    preproc.connect(reorient, 'out_file', slicer_mask2, 'in_file')
    preproc.connect(skull_strip, ('mask_file', get_overlay_args), slicer_mask2, 'overlay1')
    renamer(slicer_mask2, 'out_file', 'brain_mask_sagittal_pic')
    
    return preproc


def create_ap_freesurfer_workflow(freesurfer_dir, name="preproc_anat_freesurfer"):
    """Prepares anatomical image for fMRI analysis. This includes:
    
    1. Intensity Normalization (w/freesurfer)
    2. Skull Stripping (w/freesurfer)
    3. Reorient images (for consistency during registration)
    
    PNG images of the head and brain mask are also created for QC.
    """
    
    #####
    # Setup workflow
    #####
    
    preproc = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        "subject_id",
        "freesurfer_dir",   # $SUBJECTS_DIR
        "struct",
        "orientation"
    ]
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields), 
                            name="inputspec")
    
    # defaults
    inputnode.inputs.orientation = "RPI"    # fsl's standard brain
    
    
    #####
    # Setup output node
    #####
    
    # Outputs
    output_fields = [
        "orig",
        "head",
        "brain",
        "brain_mask",
        "head_axial_pic",
        "brain_mask_axial_pic",
        "brain_mask_sagittal_pic"        
    ]
    outputnode = pe.Node(util.IdentityInterface(fields=output_fields),
                        name="outputspec")
    
    # renaming
    renamer = SimpleOutputConnector(preproc, outputnode)
    
    
    ######
    # Commands
    ######
    
    # First call a function that may or may not run autorecon1
    skull_strip = pe.Node(run_freesurfer, name="skull_strip")
    skull_strip.inputs.directive = 'autorecon1'
    skull_strip.inputs.subjects_dir = freesurfer_dir
    preproc.connect([
        (inputnode, skull_strip, [('struct', 't1_files'),
                                  ('subject_id', 'subject_id'),
                                  ('freesurfer_dir', 'subjects_dir')])
    ])
    
    # Now get inputs
    getfree = pe.Node(nio.FreeSurferSource(subjects_dir=freesurfer_dir),
                      name="getfree")
    preproc.connect([
        (inputnode, getfree, [('freesurfer_dir', 'subjects_dir')]),
        (skull_strip, getfree, [('subject_id', 'subject_id')])
    ])
    
    # Convert to nifti
    ## orig
    convert_orig = pe.Node(fs.MRIConvert(out_type="niigz", subjects_dir=freesurfer_dir), name="convert_orig")
    preproc.connect(getfree, 'orig', convert_orig, 'in_file')
    preproc.connect([
        (inputnode, convert_orig, [('orientation', 'out_orientation'),
                                   ('freesurfer_dir', 'subjects_dir')])
    ])
    renamer(convert_orig, 'out_file', 'orig')
    ## head
    convert_head = pe.Node(fs.MRIConvert(out_type="niigz", subjects_dir=freesurfer_dir), name="convert_head")
    preproc.connect(getfree, 'T1', convert_head, 'in_file')
    preproc.connect([
        (inputnode, convert_head, [('orientation', 'out_orientation'),
                                   ('freesurfer_dir', 'subjects_dir')])
    ])
    reorient_head = pe.Node(interface=afni.Threedresample(), name='reorient_head')
    preproc.connect(convert_head, 'out_file', reorient_head, 'in_file')
    preproc.connect(inputnode, 'orientation', reorient_head, 'orientation')
    renamer(reorient_head, 'out_file', 'head')
    ## brain
    convert_brain = pe.Node(fs.MRIConvert(out_type="niigz", subjects_dir=freesurfer_dir), name="convert_brain")
    preproc.connect(getfree, 'brainmask', convert_brain, 'in_file')
    preproc.connect([
        (inputnode, convert_brain, [('orientation', 'out_orientation'),
                                    ('freesurfer_dir', 'subjects_dir')])
    ])
    reorient_brain = pe.Node(interface=afni.Threedresample(), name='reorient_brain')
    preproc.connect(convert_brain, 'out_file', reorient_brain, 'in_file')
    preproc.connect(inputnode, 'orientation', reorient_brain, 'orientation')
    renamer(reorient_brain, 'out_file', 'brain')
    
    # Create brain mask
    brain_mask = pe.Node(interface=fsl.ImageMaths(op_string='-bin'), name='brain_mask')
    preproc.connect(convert_brain, 'out_file', brain_mask, 'in_file')
    renamer(brain_mask, 'out_file', 'brain_mask')
    
    # Pic of head
    slicer_head = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                            name='slicer_head')
    preproc.connect(convert_head, 'out_file', slicer_head, 'in_file')
    renamer(slicer_head, 'out_file', 'head_axial_pic')
    
    # Pic of brain mask overlaid on head (axial)
    slicer_mask1 = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), name='slicer_mask1')
    preproc.connect(convert_head, 'out_file', slicer_mask1, 'in_file')
    preproc.connect(brain_mask, ('out_file', get_overlay_args), slicer_mask1, 'overlay1')
    renamer(slicer_mask1, 'out_file', 'brain_mask_axial_pic')
    
    # Pic of brain mask overlaid on head (sagittal)
    slicer_mask2 = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="sagittal"), name='slicer_mask2')
    preproc.connect(convert_head, 'out_file', slicer_mask2, 'in_file')
    preproc.connect(brain_mask, ('out_file', get_overlay_args), slicer_mask2, 'overlay1')
    renamer(slicer_mask2, 'out_file', 'brain_mask_sagittal_pic')
    
    return preproc


def fix_output_permissions(outputs, subject_list):
    outputs = [ op.join(outputs.basedir, s, outputs.struct) for s in subject_list ]
    for output in outputs:
        p = Process("chmod -R 775 %s" % output, to_print=True)
        if p.retcode != 0:
            print 'Error: chmod -R 775 %s' % output
    return

class AnatPreprocParser(usage.NiParser):
    def _create_parser(self, *args, **kwrds):
        """Create command-line interface"""
        parser = super(AnatPreprocParser, self)._create_parser(
            description="""
                Run preprocessing for each participant's structual images.
            """
        )
        
        group = parser.add_argument_group('Anatomical Preprocessing Options')
        group.add_argument('--struct', nargs=2, action=usage.store_io, required=True)
        group.add_argument('--run', choices=['bet', 'freesurfer'], required=True)
        group.add_argument('--orientation', default='RPI')
        group.add_argument('--freesurfer-dir', action=usage.store_directory, 
                            help="required if using --run freesurfer")
        
        # 3dinfo /Users/Shared/fsl/data/standard/MNI152_T1_2mm.nii.gz | grep orient | sed -e s/.*-orient\ // -e s/]// # to get orientation from standard brain
        return parser
    
    def _post_compile(self):
        """Check that freesurfer_dir is set if run=freesurfer"""
        super(AnatPreprocParser, self)._post_compile()
        if self.args.run == 'freesurfer' and not hasattr(self.args, 'freesurfer_dir'):
            self.parser.error("You must specify --freesurfer-dir if --run=freesurfer")
        return
    
    def _post_run(self):
        """Fix permissions of output"""
        outputs = self.args.outputs
        subject_list = self.args.subject_list
        outputs = [ op.join(outputs.basedir, s, outputs.struct) for s in subject_list ]
        for output in outputs:
            p = Process("chmod -R 775 %s" % output, to_print=True)
            if p.retcode != 0:
                print 'Error: chmod -R 775 %s' % output
        return
    

def main(arglist):
    pp = AnatPreprocParser()
    pp.run(anatomical_preprocessing, arglist)
    fix_output_permissions(pp.args.outputs, pp.args.subject_list)
    

def test_wf():
    """Testing on my computer"""
    arglist="-b /Users/zarrar/Projects/tnetworks /Users/zarrar/Projects/tnetworks/output --workingdir /Users/zarrar/Projects/tnetworks/tmp --struct *_highres.nii.gz highres --plugin Linear -s tb3417"
    main(arglist.split())

if __name__ == "__main__":
    main(sys.argv[1:])
