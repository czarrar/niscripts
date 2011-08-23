#!/usr/bin/env python
"""
This script uses nipype for generating nuisance covariates (wm, csf, global)
"""

import numpy as np

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

# my own extra stuff
import misc, usage
from utilities import SimpleOutputConnector, run_freesurfer, sink_outputs

def get_overlay_args(fname):
    """Args for the overlay1 option of slicer"""
    return (fname, 1, 1)


def regpath(regdir, fprefix):
    import os, glob
    fpath = os.path.join(regdir, fprefix)
    gpath = glob.glob(fpath)
    if gpath == 0:
        raise Exception("Could not find file '%s' in regdir '%s'" % (fprefix, regdir))
    elif gpath > 1:
        raise Exception("Too many files found for '%s' in regdir '%s'" % (fprefix, regdir))
    return gpath[0]

def fsl_prior_path(name):
    import os
    try:
        fsldir = os.environ['FSLDIR']
    except KeyError:
        raise Exception('FSL environment variables not set')
    ppath = os.path.join(fsldir, 'data', 'standard', 'tissuepriors', "avg152T1_%s.img" % name)
    if not os.path.isfile(ppath):
        raise Exception("Could not find tissue prior '%s'" % ppath)
    return ppath

tolist = lambda x: [x]


def create_tissue_masks(freesurfer_dir, name="segmentation"):
    """Generates masks of each tissue type (GM, WM, and CSF) using freesurfer.
    PNG images of each tissue type is also created for QC.
    """
    
    #####
    # Setup workflow
    #####
    
    ctissues = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        "subject_id",
        "freesurfer_dir",   # $SUBJECTS_DIR
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
        "aseg",
        "csf",
        "wm",
        "gm",
        "csf_pic",
        "wm_pic",
        "gm_pic"
    ]
    outputnode = pe.Node(util.IdentityInterface(fields=output_fields),
                        name="outputspec")
    
    # renaming
    renamer = SimpleOutputConnector(ctissues, outputnode)
    
    
    ######
    # Commands
    ######
    
    # First call a function that may or may not run autorecon1
    segment = pe.Node(run_freesurfer, name="segment")
    segment.inputs.directive = 'autorecon2'
    segment.inputs.subjects_dir = freesurfer_dir
    ctissues.connect([
        (inputnode, segment, [('subject_id', 'subject_id'),
                              ('freesurfer_dir', 'subjects_dir')])
    ])
    
    # Now get inputs
    getfree = pe.Node(nio.FreeSurferSource(subjects_dir=freesurfer_dir),
                      name="getfree")
    ctissues.connect([
        (inputnode, getfree, [('freesurfer_dir', 'subjects_dir')]),
        (segment, getfree, [('subject_id', 'subject_id')])
    ])
    
    # Save to nii.gz
    convert_aseg = pe.Node(fs.MRIConvert(out_type="niigz", subjects_dir=freesurfer_dir), 
                           name="convert_orig")
    ctissues.connect(getfree, 'aseg', convert_aseg, 'in_file')
    ctissues.connect(inputnode, 'orientation', convert_aseg, 'out_orientation')
    renamer.connect(convert_aseg, 'out_file', 'aseg')
    
    # Extract and erode a mask of the deep cerebral white matter
    extractwm = pe.Node(fs.Binarize(match=[2, 41], subjects_dir=freesurfer_dir), 
                        name="extractwm")
    ctissues.connect([
        (inputnode, extractwm, [('freesurfer_dir', 'subjects_dir')]),
        (convert_aseg, extractwm, [('out_file', 'in_file')])
    ])
    
    ## pic
    slicer_wm = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                        name='slicer_wm')
    ctissues.connect(getfree, 'T1', slicer_wm, 'in_file')
    ctissues.connect(extractwm, ('binary_file', get_overlay_args), slicer_wm, 'overlay1')
    renamer.connect(slicer_wm, 'out_file', 'wm_pic')
        
    # Extract and erode a mask of the ventricles and CSF
    extractcsf = pe.Node(fs.Binarize(match=[4, 5, 14, 15, 24, 31, 43, 44, 63], 
                                        subjects_dir=freesurfer_dir),
                         name="extractcsf")
    ctissues.connect([
        (inputnode, extractcsf, [('freesurfer_dir', 'subjects_dir')]),
        (convert_aseg, extractcsf, [('out_file', 'in_file')])
    ])
    
    ## pic
    slicer_csf = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                         name='slicer_csf')
    ctissues.connect(getfree, 'T1', slicer_csf, 'in_file')
    ctissues.connect(extractcsf, ('binary_file', get_overlay_args), slicer_csf, 'overlay1')
    renamer.connect(slicer_csf, 'out_file', 'csf_pic')
    
    # Extract a mask of the grey matter and subcortical areas and brainstem
    extractgm = pe.Node(
        fs.Binarize(match=[3,8,10,11,12,13,16,17,18,26,28,42,47,49,50,51,52,54,58,60], 
                    subjects_dir=freesurfer_dir), 
        name="extractgm"
    )
    ctissues.connect([
        (inputnode, extractgm, [('freesurfer_dir', 'subjects_dir')]),
        (convert_aseg, extractgm, [('out_file', 'in_file')])
    ])
    
    ## pic
    slicer_gm = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                         name='slicer_gm')
    ctissues.connect(getfree, 'T1', slicer_gm, 'in_file')
    ctissues.connect(extractgm, ('binary_file', get_overlay_args), slicer_gm, 'overlay1')
    renamer.connect(slicer_gm, 'out_file', 'gm_pic')
    
    return ctissues


def create_nuisance_mask_workflow(fwhm, mask_fprefix, freesurfer_dir, name="nuisance_mask"):
    """Transforms a csf/wm/gm mask from anatomical to functional space and then some.
    """
    
    #####
    # Setup workflow
    #####
    
    nuisance = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        "mask_file",
        "prior_file",
        "reg_dir",
        "brain_mask_file",
        "erode", 
        "threshold",
        "freesurfer_dir"
    ]
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields), 
                            name="inputspec")
    
    inputnode.inputs.threshold = 0.75
    
    
    #####
    # Setup output node
    #####
    
    # Outputs
    output_fields = [
        "mask",
        "mask_pic"
    ]
    outputnode = pe.Node(util.IdentityInterface(fields=output_fields),
                        name="outputspec")
    
    # renaming
    renamer = SimpleOutputConnector(nuisance, outputnode)
    
    
    #####
    # Commands
    #####
    
    # 1. erode
    erode = pe.Node(fs.Binarize(match=[1], subjects_dir=freesurfer_dir), name="01_erode")
    nuisance.connect([
        (inputnode, erode, [('mask_file', 'in_file'),
                            ('erode', 'erode'),
                            ('freesurfer_dir', 'subjects_dir')])
    ])
    
    # 2. anat space => func space
    mask2func = pe.Node(fsl.ApplyXfm(apply_xfm=True, interp="trilinear"), 
                        name="02_mask2func")
    nuisance.connect([
        (erode, mask2func, [('out_file', 'in_file')]),
        (inputnode, mask2func, [(('reg_dir', regpath, 'func.*'), 'reference'),
                                (('reg_dir', regpath, 'highres2func.mat'), 'in_matrix_file')])
    ])
    
    # 3. smooth
    if fwhm > 0:
        smooth = pe.Node(fsl.IsotropicSmooth(fwhm=fwhm), 
                         name="03_smooth")
        nuisance.connect([
            (mask2func, smooth, [('out_file', 'in_file')])
        ])
        nextnode = smooth
    else:
        nextnode = anat2func
    
    # 4. prior: std space => func space
    prior2func = pe.Node(fsl.ApplyXfm(apply_xfm=True, interp="nearestneighbour"), 
                     name="04_prior2func")
    nuisance.connect([
        (inputnode, prior2func, [('prior_file', 'in_file'), 
                                 (('reg_dir', regpath, 'func.*'), 'reference'),
                                 (('reg_dir', regpath, 'standard2func.mat'), 'in_matrix_file')])
    ])
    
    # 5. mask by prior
    prior_masked = pe.Node(fsl.ApplyMask(),
                           name="05_prior_masked")
    nuisance.connect([
        (nextnode, prior_masked, [('out_file', 'in_file')]),
        (prior2func, prior_masked, [('out_file', 'mask_file')])
    ])
    
    # 6. threshold
    threshold = pe.Node(fsl.Threshold(direction="below"),
                        name="06_threshold")
    nuisance.connect(inputnode, 'threshold', threshold, 'thresh')
    nuisance.connect(prior_masked, 'out_file', threshold, 'in_file')
    
    # 7. binarize
    binarize = pe.Node(fsl.UnaryMaths(operation="bin"),
                        name="07_binarize")
    nuisance.connect(threshold, 'out_file', binarize, 'in_file')
    
    # 8. mask by brain
    brain_masked = pe.Node(fsl.ApplyMask(),
                           name="08_brain_masked")
    nuisance.connect([
        (binarize, brain_masked, [('out_file', 'in_file')]),
        (inputnode, brain_masked, [('brain_mask_file', 'mask_file')])
    ])
    renamer.connect(brain_masked, 'out_file', 'mask', format_string=mask_fprefix)
    
    # 9. picture of mask over func
    slicer_mask = pe.Node(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                         name='09_slicer_mask')
    nuisance.connect(inputnode, ('reg_dir', regpath, 'func.*'), slicer_mask, 'in_file')
    nuisance.connect(brain_masked, ('out_file', get_overlay_args), slicer_mask, 'overlay1')
    renamer.connect(slicer_mask, 'out_file', 'mask_pic', format_string=mask_fprefix)
    
    return nuisance


def create_nuisance_evs_workflow(freesurfer_dir, fwhm, name="nuisance_evs"):
    
    #####
    # Setup workflow
    #####
    
    wf = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        # segment
        "subject_id",       # required
        "freesurfer_dir", 
        "orientation", 
        # mask creation
        "brain_mask",       # required
        "reg_dir",          # required
        "csf_prior", 
        "csf_erode", 
        "csf_threshold", 
        "wm_prior", 
        "wm_erode", 
        "wm_threshold",
        # ts extraction 
        "func",             # required
    ]
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields), 
                            name="inputspec")
    
    # defaults
    inputnode.inputs.freesurfer_dir = freesurfer_dir
    inputnode.inputs.orientation = 'RPI'
    inputnode.inputs.csf_prior = fsl_prior_path("csf")
    inputnode.inputs.csf_erode = 0
    inputnode.inputs.csf_threshold = 0.75
    inputnode.inputs.wm_prior = fsl_prior_path("white")
    inputnode.inputs.wm_erode = 2
    inputnode.inputs.wm_threshold = 0.75
        
    
    #####
    # Setup output node
    #####
    
    # Outputs
    output_fields_highres = [
        # segmentation
        "csf",
        "wm",
        "gm",
        "csf_pic",
        "wm_pic",
        "gm_pic"
    ]
    outputnode_highres = pe.Node(util.IdentityInterface(fields=output_fields_highres),
                                                        name="outputspec_highres")
    
    output_fields_func = [
        # masks
        "mask_csf",
        "mask_wm",
        "mask_csf_pic", 
        "mask_wm_pic", 
        # ts
        "ts_global",
        "ts_csf",
        "ts_wm",
        "nuisance_tsplot_pic"
    ]
    outputnode_func = pe.Node(util.IdentityInterface(fields=output_fields_func),
                                                     name="outputspec_func")
    
    # renaming
    renamer_highres = SimpleOutputConnector(wf, outputnode_highres)
    renamer_func = SimpleOutputConnector(wf, outputnode_func)
    
    
    #####
    # COMMANDS
    #####
    
    # Segment Brain
    segment = create_tissue_masks(freesurfer_dir, name="01_segment")
    wf.connect([
        (inputnode, segment, [('subject_id', 'inputspec.subject_id'),
                              ('freesurfer_dir', 'inputspec.freesurfer_dir'),
                              ('orientation', 'inputspec.orientation')]),
        (segment, outputnode_highres, [('outputspec.csf', 'csf'),
                                       ('outputspec.wm', 'wm'),
                                       ('outputspec.gm', 'gm'), 
                                       ('outputspec.csf_pic', 'csf_pic'), 
                                       ('outputspec.wm_pic', 'wm_pic'),
                                       ('outputspec.gm_pic', 'gm_pic')])
    ])
    
    # Create CSF mask
    csfmask = create_nuisance_mask_workflow(name="02_mask_csf", fwhm=fwhm, mask_fprefix='mask_csf', 
                                                freesurfer_dir=freesurfer_dir)
    wf.connect([
        (segment, csfmask, [('outputspec.csf', 'inputspec.mask_file')]), 
        (inputnode, csfmask, [('brain_mask', 'inputspec.brain_mask_file'),
                              ('csf_prior', 'inputspec.prior_file'),
                              ('reg_dir', 'inputspec.reg_dir'),
                              ('freesurfer_dir', 'inputspec.freesurfer_dir')]),
        (csfmask, outputnode_func, [('outputspec.mask', 'mask_csf'),
                                    ('outputspec.mask_pic', 'mask_csf_pic')])
    ])
    
    # Create WM mask
    wmmask = create_nuisance_mask_workflow(name="02_mask_wm", fwhm=fwhm, mask_fprefix='mask_wm', 
                                            freesurfer_dir=freesurfer_dir)
    wf.connect([
        (segment, wmmask, [('outputspec.wm', 'inputspec.mask_file')]), 
        (inputnode, wmmask, [('brain_mask', 'inputspec.brain_mask_file'),
                              ('wm_prior', 'inputspec.prior_file'),
                              ('reg_dir', 'inputspec.reg_dir'),
                              ('freesurfer_dir', 'inputspec.freesurfer_dir')]),
        (wmmask, outputnode_func, [('outputspec.mask', 'mask_wm'),
                                   ('outputspec.mask_pic', 'mask_wm_pic')])
    ])
    
    # Extract TS for global
    meants_global = pe.Node(fsl.ImageMeants(), 
                            name="03_meants_global")
    wf.connect([
        (inputnode, meants_global, [('func', 'in_file'),
                                     ('brain_mask', 'mask_file')])
    ])
    renamer_func.connect(meants_global, 'out_file', 'ts_global')
    renamer_func.connect(inputnode, 'brain_mask', 'mask_global')
    
    # Extract TS for csf
    meants_csf = pe.Node(fsl.ImageMeants(), 
                            name="03_meants_csf")
    wf.connect([
        (inputnode, meants_csf, [('func', 'in_file')]),
        (csfmask, meants_csf, [('mask', 'mask_file')]),
    ])
    renamer_func.connect(meants_csf, 'out_file', 'ts_csf')
    
    # Extract TS for wm
    meants_csf = pe.Node(fsl.ImageMeants(), 
                            name="03_meants_wm")
    wf.connect([
        (inputnode, meants_wm, [('func', 'in_file')]),
        (wmmask, meants_wm, [('mask', 'mask_file')]),
    ])
    renamer_func.connect(meants_wm, 'out_file', 'ts_wm')
    
    # Create plot
    ## concat
    concatnode = pe.Node(interface=util.Merge(3),
                         name='04_concat')
    wf.connect(meants_global, ('out_file', tolist), concatnode, 'in1')
    wf.connect(meants_csf, ('out_file', tolist), concatnode, 'in2')
    wf.connect(meants_wm, ('out_file', tolist), concatnode, 'in3')
    ## plot
    tsplot = pe.Node(interface=fsl.PlotTimeSeries(title='Nuisance Time-Series', 
                                                  labels=['global', 'csf', 'wm']), 
                     name="04_tsplot")
    wf.connect(concatnode, 'out', tsplot, 'in_file')
    renamer_func.connect(tsplot, 'out_file', 'nuisance_tsplot_pic')
    
    return wf


def nuisance_evs(
    subject_list, 
    inputs, outputs, workingdir, output_type, 
    orientation, fwhm, 
    name="nuisance_evs"):
    """Wrapper...
    """
    
    #####
    # Setup workflow
    #####
    
    wf = create_nuisance_evs_workflow(inputs.freesurfer_dir, fwhm, name=name)
    wf.base_dir = workingdir
    
    ######
    # Setup data source
    ######
    
    # Subject source node
    subinfo = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name='subinfo', 
                        iterables=('subject_id', subject_list))
    
    # Say where to find input data
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['func', 'func_mask']), name='datasource')
    datasource.inputs.base_directory=os.path.abspath(inputs.basedir)
    datasource.inputs.template = "*"
    datasource.inputs.field_template = dict(
        func = os.path.join("%s", inputs.func),
        func_mask = os.path.join("%s", inputs.func_mask)
    )
    datasource.inputs.template_args = dict(
                                        func = [['subject_id']],
                                        func_mask = [['subject_id']]
                                      )
    
    # Link inputs
    wf.connect([
        (subinfo, datasource, [('subject_id', 'subject_id')]), 
        (subinfo, wf, [('subject_id', 'inputspec.subject_id')]), 
        (datasource, wf, [('func', 'inputspec.func'),
                          ('func_mask', 'inputspec.brain_mask')])
    ])
    wf.inputs.inputspec.reg_dir = inputs.reg_dir
    
    
    ######
    # Setup data sink
    ######
    
    datasink_highres = pe.Node(interface=nio.DataSink(), name='datasink_highres')
    datasink_highres.inputs.base_directory = os.path.abspath(outputs.basedir)
    datasink_highres.inputs.regexp_substitutions = (r"_subject_id_(\w|\d)+", outputs.segment)
    wf.connect(subinfo, 'subject_id', datasink_highres, 'container')
    sink_outputs(wf, wf.get_node("outputspec_highres"), datasink_highres)
    
    datasink_func = pe.Node(interface=nio.DataSink(), name='datasink_func')
    datasink_func.inputs.base_directory = os.path.abspath(outputs.basedir)
    datasink_func.inputs.regexp_substitutions = (r"_subject_id_(\w|\d)+", outputs.nuisance)
    wf.connect(subinfo, 'subject_id', datasink_func, 'container')
    sink_outputs(wf, wf.get_node("outputspec_func"), datasink_func)
    
    return wf


class NuisanceEVsParser(usage.NiParser):
    def _create_parser(self, *args, **kwrds):
        """Create command-line interface"""
        parser = super(NuisanceEVsParser, self)._create_parser(
            description="""
                Create nuisance covariates (wm, csf, and global signal)
            """
        )
        
        group = parser.add_argument_group('Nuisance EV Creation Options')
        group.add_argument('--func', action=usage.store_input, required=True)
        group.add_argument('--func-mask', action=usage.store_input, required=True)
        group.add_argument('--reg-dir', action=usage.store_input, required=True)
        group.add_argument('--freesurfer-dir', action=usage.store_input, check_dir=True, 
                           required=True)
        group.add_argument('--segment', action=usage.store_output, required=True)
        group.add_argument('--nuisance', action=usage.store_output, required=True)
        group.add_argument('--fwhm', type=float, required=True)
        group.add_argument('--orientation', default='RPI')
        
        # todo: can auto determine orientation desired based on standard space to registered to
        # 3dinfo /Users/Shared/fsl/data/standard/MNI152_T1_2mm.nii.gz | grep orient | sed -e s/.*-orient\ // -e s/]// # to get orientation from standard brain
        
        return parser
    
    def _post_run(self):
        """Fix permissions of output"""
        outputs = self.args.outputs
        subject_list = self.args.subject_list
        outputs1 = [ op.join(outputs.basedir, s, outputs.segment) for s in subject_list ]
        outputs2 = [ op.join(outputs.basedir, s, outputs.nuisance) for s in subject_list ]
        outputs = outputs1 + outputs2
        for output in outputs:
            p = Process("chmod -R 775 %s" % output, to_print=True)
            if p.retcode != 0:
                print 'Error: chmod -R 775 %s' % output
        return
    

def main(arglist):
    pp = NuisanceEVsParser()
    pp.run(nuisance_evs, arglist)


def test_wf():
    """Testing on my computer"""
    arglist="-b /Users/zarrar/Projects/tnetworks /Users/zarrar/Projects/tnetworks/output --workingdir /Users/zarrar/Projects/tnetworks/tmp --struct *_highres.nii.gz highres --plugin Linear -s tb3417"
    main(arglist.split())

if __name__ == "__main__":
    main(sys.argv[1:])
