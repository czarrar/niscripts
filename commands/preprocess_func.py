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
import nipype.algorithms.rapidart as ra
import nipype.interfaces.utility as util # utility
import nipype.pipeline.engine as pe # pypeline engine
import os # system functions
import argparse # command-line
import os.path as op
from nipype.interfaces.base import traits

import sys
sys.path.append(op.join(os.environ.get("NISCRIPTS"), "include"))

import misc, usage # my own extra stuff
from utilities import *
from execute import Process

from nipype.workflows.fsl import create_susan_smooth



def pickfirst(files):
    if isinstance(files, list):
        return files[0]
    else:
        return files

def pickmiddle(files):
    from nibabel import load
    import numpy as np
    middlevol = []
    for f in files:
        middlevol.append(int(np.ceil(load(f).get_shape()[3]/2)))
    return middlevol

def pickvol(filenames, fileidx, which):
    from nibabel import load
    import numpy as np
    if which.lower() == 'first':
        idx = 0
    elif which.lower() == 'middle':
        idx = int(np.ceil(load(filenames[fileidx]).get_shape()[3]/2))
    else:
        raise Exception('unknown value for volume selection : %s'%which)
    return idx

def pick_list_elem(l, index):
    return l[index]

# borrowed from Fluid_NiPype
def max_motion_func(rms_files):
    """Determine the maximum absolute and relative motion values."""
    from os import getcwd
    from os.path import join
    from numpy import loadtxt, max
    motion = map(loadtxt, rms_files)
    maxima = map(max, motion)
    out_file = join(getcwd(), "max_motion.txt")
    with open(out_file, "w") as f:
        f.write("#Absolute:\n%.4f\n#Relative\n%.4f"%tuple(maxima))
    return out_file

def get_thresh_op(thresh):
    return ['-thr %.10f -Tmin -bin'%(0.1*val[1]) for val in thresh]

def chooseindex(fwhm):
    if fwhm<1:
        return [0]
    else:
        return [1]

def fwhm_used(fwhm):
    if fwhm<1:
        return str(0)
    else:
        return str(int(fwhm))

def filter_used(filt):
    if filt<0:
        return '0'
    else:
        return '1'

def getmeanscale(medianvals):
    return ['-mul %.10f'%(10000./val) for val in medianvals]

def combine_motion_params(in_files):
    import sys
    from os import getcwd, environ
    from os.path import join, basename, splitext
    sys.path.append(join(environ.get("NISCRIPTS"), "include"))
    from execute import Process
    
    cwd = getcwd()
    ext = "png"
    out_file = join(cwd, basename(in_files[0]))
    f = file(out_file, 'w')
    
    tmp = Process("cat %s" % " ".join(in_files), stdout=f, cwd=cwd, to_print=True)
    print tmp.stderr
    
    f.close()
    return out_file

def split_motion_params(in_file, out_prefix="motion"):
    import numpy as np
    import os
    import os.path as op
    
    cwd = os.getcwd()
    xmat = np.loadtxt(in_file)
    nc = xmat.shape[1]
    out_files = [ op.join(cwd, "%s_%02i.1D" % (out_prefix, i+1)) for i in xrange(nc) ]
    for i in xrange(nc):
        xcol = [ str(x) for x in xmat[:,i].tolist() ]
        f = file(out_files[i], 'w')
        f.write("\n".join(xcol))
        f.close()
    
    return out_files

def get_overlay_args(fname):
    """Args for the overlay1 option of slicer"""
    return (fname, 1, 1)

tolist = lambda x: [x]



def functional_preprocessing(
    subject_list, 
    inputs, outputs, workingdir, output_type, 
    fwhm, hpfilter, lpfilter, tr, motion_nstages, label, 
    name="%s_preprocessing", whichvol="middle", timeshift=False, tpattern=None):
    """todo"""
    
    #####
    # Setup pipeline
    #####
    
    if name.find("%s") != -1:
        name = name % label
    else:
        print "ERROR: You must have a '%s' in the name for the label"
        raise SystemExit(2)
    
    preproc = create_func_preproc_workflow(name, whichvol, timeshift, tpattern)
    preproc.base_dir = workingdir
    
    preproc.inputs.inputspec.fwhm = fwhm    
    preproc.inputs.inputspec.highpass = float(hpfilter)/(2*tr)
    preproc.inputs.inputspec.lowpass = float(lpfilter)/(2*tr)
    preproc.inputs.inputspec.motion_nstages = motion_nstages
    
    ## save for later
    inputnode = preproc.get_node("inputspec")
    outputnode = preproc.get_node("outputspec")
    
    
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
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], 
                                                    outfields=['func'],
                                                    sort_filelist=True), 
                         name='datasource')
    datasource.inputs.base_directory=os.path.abspath(inputs.basedir)
    datasource.inputs.template = os.path.join("%s", inputs.func)
    datasource.inputs.template_args = dict(func=[['subject_id']])
    
    # Get the number of runs for each participant (for renaming purposes)
    datasource.inputs.subject_id = subject_list
    ds = datasource.run()
    runs = [ len(x) for x in ds.outputs.func ]  # list with number of runs per subject
    max_runs = max(runs)
    
    # Link inputs
    preproc.connect(subinfo, 'subject_id', datasource, 'subject_id')
    preproc.connect(datasource, 'func', inputnode, 'func')
    
    
    ######
    # Setup data sink
    ######
    
    # Datasink
    ## will get: "base_directory/subject_id/output_anatdir"
    datasink = pe.Node(interface=nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = outputs.basedir
    
    ## substitute map stuff
    preprocsinksubs = get_mapnode_substitutions(preproc, outputnode, max_runs, 
        unmap=["func_mc_ref", "example_func_all", "func_mask_all", "motion_all", 
               "motion_01", "motion_02", "motion_03", "motion_04", "motion_05", "motion_06"])
    datasink.inputs.substitutions = preprocsinksubs
    
    # replace subject_id stuff with functional scan
    datasink.inputs.regexp_substitutions = (r"_subject_id_(\w|\d)+", outputs.func)
    
    ## connect
    preproc.connect(subinfo, 'subject_id', datasink, 'container')
    output_fields = outputnode.outputs.get()
    for field in output_fields:
        preproc.connect(outputnode, field, datasink, "@%s" % field)
    
    return preproc

# whichvol can be first, middle, or mean
# hp should be in TRs
def create_func_preproc_workflow(name='functional_preprocessing', whichvol='middle', timeshift=False, tpattern=None):
        
    #####
    # Setup workflow
    #####
    
    preproc = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        "func",     # list of functional images to run
        "fwhm",     # fwhm
        "highpass", # highpass filtering
        "lowpass",  # lowpass filtering
        "motion_nstages"
    ]
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields), 
                            name="inputspec")
    
    
    #####
    # Setup output node
    #####
    
    # Outputs
    output_fields = [
        "func_mc_ref",      # reference for motion correction
        "motion",           # motion parameters
        "motion_rot",       # plot of rotation
        "motion_trans",     # plot of translations
        "motion_disp",      # plot of abs/rel displacement
        "motion_max",       # textfile with max abs/rel displacement
        "func_mc",          # 4D motion corrected and skull stripped data
        "func_afnimask",    # mask from 3dAutomask
        "example_func",     # useful for registration later (simply a mean image)
        "func_mask",        # mask that is less contrained
        "func_mc_sm",       # smoothed data
        "func_mc_sm_ft",    # time filtered data
        "func_preproc",     # final output (also has been intensity normalized)
        "func_mean",        # mean of final output
        "example_func_all", 
        "func_mask_all",
        "motion_all",
        "motion_01",
        "motion_02",
        "motion_03",
        "motion_04",
        "motion_05",
        "motion_06",
        "pics_func_mean_head1",
        "pics_func_mean_head2",
        "pics_func_mean_brain1",
        "pics_func_mean_brain2"
    ]
    outputnode = pe.Node(util.IdentityInterface(fields=output_fields),
                        name="outputspec")
    
    #####
    # Setup renaming
    #####
    
    renamer = OutputConnector(preproc, outputnode)  # allows easy renaming of output file names
    
    
    #####
    # Main Commands
    ####
    
    """
    Prepare Input
    """
    
    # Convert functional images to float representation
    img2float = pe.MapNode(fsl.ChangeDataType(output_datatype="float"),
                           iterfield=["in_file"],
                           name="00_img2float")
    preproc.connect(inputnode, 'func', img2float, 'in_file')
    nextnode = img2float
    
    """
    Time Shift Slices
    """
    
    if timeshift:
        tshift = pe.MapNode(afni.ThreedTshift(tzero=0),
                            iterfield=["in_file"], name="01_tshift")
        if tpattern is not None:
            tshift.inputs.tpattern = tpattern
        preproc.connect(img2float, 'out_file', tshift, 'in_file')
        nextnode = tshift
    
    """
    Deoblique and Reorient to FSL Friendly Space
    """
    
    deoblique = pe.MapNode(interface=afni.ThreedWarp(deoblique=True), 
                            iterfield=["in_file"], name='01_deoblique')
    preproc.connect(nextnode, "out_file", deoblique, "in_file")
    
    # TODO:
    # get orientation using following command:
    # orient=$( 3dinfo /Users/Shared/fsl/data/standard/MNI152_T1_2mm.nii.gz | grep orient | sed -e s/.*orient// -e s/\]// )
    # so add additional input of reference orientation! (maybe this can be a general option?)
    reorient = pe.MapNode(interface=afni.Threedresample(orientation='RPI'), 
                            iterfield=["in_file"], name='01_reorient')
    preproc.connect(deoblique, "out_file", reorient, "in_file")
    
    """
    Motion Correction
    """
    
    # Get the middle volume of each run for motion correction 
    if whichvol != 'mean':
        extract_ref = pe.Node(fsl.ExtractROI(t_size=1), iterfield=["in_file", "t_min"],
                                name = "02_extractref")
        preproc.connect(reorient, ('out_file', pickfirst), extract_ref, 'in_file')
        preproc.connect(reorient, ('out_file', pickvol, 0, whichvol), extract_ref, 't_min')
        renamer.connect(extract_ref, 'roi_file', 'func_mc_ref')
    
    # Realign the functional runs to the reference (some volume from the 1st run)
    motion_correct = pe.MapNode(interface=fsl.MCFLIRT(save_mats = True,
                                                      save_plots = True,
                                                      save_rms = True,
                                                      interpolation = 'sinc'),
                                name='02_realign', iterfield = ['in_file'])
    preproc.connect(reorient, 'out_file', motion_correct, 'in_file')
    preproc.connect(inputnode, 'motion_nstages', motion_correct, 'stages')
    if whichvol != 'mean':
        preproc.connect(extract_ref, 'roi_file', motion_correct, 'ref_file')
    else:
        motion_correct.inputs.mean_vol = True
        renamer.connect(motion_correct, ('mean_img', pickfirst), 'func_mc_ref')
    renamer.connect(motion_correct, 'par_file', 'motion')
    
    # Combine motion parameters from different runs
    combinemotion = pe.Node(util.Function(input_names=["in_files"],
                                    output_names=["out_file"],
                                    function=combine_motion_params),
                           name="combinemotion")
    preproc.connect(motion_correct, 'par_file', combinemotion, 'in_files')
    renamer.connect(combinemotion, 'out_file', 'motion_all')
    
    # Split motion up
    splitmotion = pe.Node(util.Function(input_names=["in_file", "out_prefix"],
                                        output_names=["out_files"], 
                                        function=split_motion_params), 
                          name="splitmotion")
    preproc.connect(combinemotion, "out_file", splitmotion, "in_file")
    preproc.connect(splitmotion, ("out_files", pick_list_elem, 0), outputnode, 'motion_01')
    preproc.connect(splitmotion, ("out_files", pick_list_elem, 1), outputnode, 'motion_02')
    preproc.connect(splitmotion, ("out_files", pick_list_elem, 2), outputnode, 'motion_03')
    preproc.connect(splitmotion, ("out_files", pick_list_elem, 3), outputnode, 'motion_04')
    preproc.connect(splitmotion, ("out_files", pick_list_elem, 4), outputnode, 'motion_05')
    preproc.connect(splitmotion, ("out_files", pick_list_elem, 5), outputnode, 'motion_06')
    
    # Plot rotation parameters from MCFLIRT
    plotrot = pe.MapNode(fsl.PlotMotionParams(in_source="fsl", plot_type="rotations"),
                         name="02_plotrotation", iterfield=["in_file"])
    preproc.connect(motion_correct, 'par_file', plotrot, 'in_file')
    renamer.connect(plotrot, 'out_file', 'motion_rot')
    
    # Plot translation parameters from MCFLIRT
    plottrans = pe.MapNode(fsl.PlotMotionParams(in_source="fsl", plot_type="translations"),
                           name="02_plottranslation", iterfield=["in_file"])
    preproc.connect(motion_correct, 'par_file', plottrans, 'in_file')
    renamer.connect(plottrans, 'out_file', 'motion_trans')
    
    # Plot displacement parameters from MCFLIRT
    plotdisp = pe.MapNode(fsl.PlotMotionParams(in_source="fsl", plot_type="displacement"),
                          name="02_plotdisplacement", iterfield=["in_file"])
    preproc.connect(motion_correct, 'rms_files', plotdisp, 'in_file')
    renamer.connect(plotdisp, 'out_file', 'motion_disp')
    
    # Plot maximum displacement parameters (abs/rel) from MCFLIRT
    maxmotion = pe.MapNode(util.Function(input_names=["rms_files"],
                                         output_names=["out_file"],
                                         function=max_motion_func),
                           iterfield=["rms_files"], name="02_maxmotion")
    preproc.connect(motion_correct, 'rms_files', maxmotion, 'rms_files')
    renamer.connect(maxmotion, 'out_file', 'motion_max')
    
    """
    Skull Strip
    """
    
    # Get a mean image of the realigned timeseries
    meanfunc1 = pe.MapNode(fsl.MeanImage(),
                           iterfield=["in_file"],
                           name="03_meanfunc1")
    preproc.connect(motion_correct, 'out_file', meanfunc1, 'in_file')
    
    # Get slices
    slicer_head1 = pe.MapNode(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                                iterfield = ["in_file"], 
                                name = 'slicer_head1')
    slicer_head2 = pe.MapNode(interface=misc.Slicer(width=5, height=4, slice_name="sagittal"), 
                                iterfield = ["in_file"], 
                                name = 'slicer_head2')
    preproc.connect([
        (meanfunc1, slicer_head1, [('out_file', 'in_file')]),
        (meanfunc1, slicer_head2, [('out_file', 'in_file')])
    ])
    renamer.connect(slicer_head1, "out_file", "pics_func_mean_head1")
    renamer.connect(slicer_head2, "out_file", "pics_func_mean_head2")
    
    # Skullstrip the mean functional image
    meanmask1 = pe.MapNode(afni.ThreedAutomask(dilate = 1),
                           iterfield = ["in_file"],
                           name = "03_meanmask1")
    preproc.connect(meanfunc1, 'out_file', meanmask1, 'in_file')
    renamer.connect(meanmask1, 'out_file', 'func_afnimask')
    
    # Apply to mean
    meanfunc2 = pe.MapNode(fsl.ApplyMask(),
                            iterfield=["in_file", "mask_file"],
                            name = "03_meanfunc2")
    preproc.connect(meanfunc1, 'out_file', meanfunc2, 'in_file')
    preproc.connect(meanmask1, 'out_file', meanfunc2, 'mask_file')
    renamer.connect(meanfunc2, 'out_file', 'example_func')
    
    # Get slices
    slicer_brain1 = pe.MapNode(interface=misc.Slicer(width=5, height=4, slice_name="axial"), 
                                iterfield = ["in_file"], 
                                name = 'slicer_brain1')
    slicer_brain2 = pe.MapNode(interface=misc.Slicer(width=5, height=4, slice_name="sagittal"), 
                                iterfield = ["in_file"], 
                                name = 'slicer_brain2')
    preproc.connect([
        (meanfunc2, slicer_brain1, [('out_file', 'in_file')]),
        (meanfunc2, slicer_brain2, [('out_file', 'in_file')])
    ])
    renamer.connect(slicer_brain1, "out_file", "pics_func_mean_brain1")
    renamer.connect(slicer_brain2, "out_file", "pics_func_mean_brain2")
    
    # Combine different means and get their mean
    mergenode = pe.Node(fsl.Merge(dimension="t"), name="03b_merge")
    preproc.connect(meanfunc2, 'out_file', mergenode, 'in_files')
    meanfunc = pe.Node(fsl.MeanImage(), name="03b_meanfunc")
    preproc.connect(mergenode, 'merged_file', meanfunc, 'in_file')
    renamer.connect(meanfunc, 'out_file', 'example_func_all', format_string="example_func")
    
    # Apply to 4D functional
    funcbrain1 = pe.MapNode(fsl.ApplyMask(),
                           iterfield=["in_file", "mask_file"],
                           name = "03_funcbrain1")
    preproc.connect(motion_correct, 'out_file', funcbrain1, 'in_file')
    preproc.connect(meanmask1, 'out_file', funcbrain1, 'mask_file')
    renamer.connect(funcbrain1, 'out_file', 'func_mc')
    
    """
    Dilated Brain Mask
    """
    
    # Determine the 2nd and 98th percentile intensities of each run
    getthresh = pe.MapNode(fsl.ImageStats(op_string="-p 2 -p 98"),
                           iterfield = ["in_file"],
                           name="04_getthreshold")
    preproc.connect(funcbrain1, 'out_file', getthresh, 'in_file')
    
    # Threshold the functional data at 10% of the 98th percentile
    threshold = pe.MapNode(fsl.ImageMaths(out_data_type="char",
                                          suffix="_thresh"),
                           iterfield = ["in_file", "op_string"],
                           name="04_threshold")
    preproc.connect(funcbrain1, 'out_file', threshold, "in_file")
    preproc.connect(getthresh, ("out_stat", get_thresh_op), threshold, "op_string")
    
    # Determine the median value of the functional runs using the mask
    medianval = pe.MapNode(fsl.ImageStats(op_string="-k %s -p 50"),
                           iterfield = ["in_file", "mask_file"],
                           name="04_medianval")
    preproc.connect(motion_correct, 'out_file', medianval, 'in_file')
    preproc.connect(threshold, 'out_file', medianval, 'mask_file')
    
    # Dilate the mask
    dilatemask = pe.MapNode(fsl.DilateImage(operation="max"),
                            iterfield=["in_file"],
                            name="04_dilatemask")
    preproc.connect(threshold, 'out_file', dilatemask, 'in_file')
    renamer.connect(dilatemask, 'out_file', 'func_mask')
    
    # Combine masks
    mergenode = pe.Node(fsl.Merge(dimension="t"), name="04b_merge")
    preproc.connect(dilatemask, 'out_file', mergenode, 'in_files')
    maskfunc = pe.Node(fsl.ImageMaths(op_string="-Tmin", 
                                          suffix="_mask"),
                           name="04b_maskfunc")
    preproc.connect(mergenode, 'merged_file', maskfunc, 'in_file')
    renamer.connect(maskfunc, 'out_file', 'func_mask_all', format_string="func_mask")
    
    # Mask the runs again with this new mask
    funcbrain2 = pe.MapNode(fsl.ApplyMask(),
                           iterfield=["in_file", "mask_file"],
                           name="04_funcbrain2")
    preproc.connect(motion_correct, 'out_file', funcbrain2, 'in_file')
    preproc.connect(dilatemask, 'out_file', funcbrain2, 'mask_file')
    
    # Get a new mean image from each functional run
    meanfunc3 = pe.MapNode(fsl.MeanImage(),
                           iterfield=["in_file"],
                           name="04_meanfunc3")
    preproc.connect(funcbrain2, 'out_file', meanfunc3, 'in_file')        
    
    """
    Detect Any Artifacts
    """
    
    art = pe.MapNode(ra.ArtifactDetect(use_differences = [True, False],
                                       use_norm = True,
                                       zintensity_threshold = 3,
                                       norm_threshold = 1,
                                       parameter_source = "FSL",
                                       mask_type = "file"),
                     iterfield=["realignment_parameters","realigned_files","mask_file"],
                     name="05_art")
    preproc.connect(motion_correct, "par_file", art, "realignment_parameters")
    preproc.connect(funcbrain2, "out_file", art, "realigned_files")
    preproc.connect(dilatemask, "out_file", art, "mask_file")
    
    plotmean = pe.MapNode(fsl.PlotTimeSeries(title="Global Mean Intensity"),
                          iterfield=["in_file"],
                          name="05_plotmean")
    preproc.connect(art, 'intensity_files', plotmean, 'in_file')
    
    """
    Smoothing
    """
    
    # Anisitropic-like? smoothing with your friend susan
    smooth = create_susan_smooth(name="06_susan_smooth")
    preproc.connect(inputnode, 'fwhm', smooth, 'inputnode.fwhm')
    preproc.connect(funcbrain2, 'out_file', smooth, 'inputnode.in_files')
    preproc.connect(dilatemask, 'out_file', smooth, 'inputnode.mask_file')
    
    # Mask smoothed data with dilated mask
    funcbrain3 = pe.MapNode(fsl.ApplyMask(),
                           iterfield=["in_file", "mask_file"],
                           name="06_funcbrain3")
    preproc.connect(smooth, 'outputnode.smoothed_files', funcbrain3, 'in_file')
    preproc.connect(dilatemask, 'out_file', funcbrain3, 'mask_file')
    
    # Determine if want to take output from smoothed or non-smoothed data forward
    # This would be the case if the fwhm is less than 1/3 of the voxel size
    ## gather 2 types of functional images
    concatnode = pe.Node(interface=util.Merge(2),
                         name='06_concat')
    preproc.connect(funcbrain2, ('out_file', tolist), concatnode, 'in1')
    preproc.connect(funcbrain3, ('out_file', tolist), concatnode, 'in2')
    ## select one
    selectnode = pe.Node(interface=util.Select(),name='06_select')
    preproc.connect(concatnode, 'out', selectnode, 'inlist')
    preproc.connect(inputnode, ('fwhm', chooseindex), selectnode, 'index')
    rename_smooth = renamer.connect(selectnode, 'out', 'func_mc_sm', 
                                    format_string="func_mc_sm%(fwhm)s", 
                                    fwhm=None)
    preproc.connect(inputnode, ('fwhm', fwhm_used), rename_smooth, 'fwhm')
    
    """
    Filter
    """
    
    filt = pe.MapNode(fsl.TemporalFilter(),
                          iterfield=["in_file"],
                          name="07_filter")
    preproc.connect(inputnode, 'highpass', filt, 'highpass_sigma')
    preproc.connect(inputnode, 'lowpass', filt, 'lowpass_sigma')
    preproc.connect(selectnode, 'out', filt, 'in_file')
    ## set renamed output
    rename_filt = renamer.connect(filt, 'out_file', 'func_mc_sm_ft', 
                                  format_string="func_mc_sm%(fwhm)s_hp%(hp)s_lp%(lp)s",
                                  fwhm=None, hp=None, lp=None)
    preproc.connect(inputnode, ('fwhm', fwhm_used), rename_filt, 'fwhm')
    preproc.connect(inputnode, ('highpass', filter_used), rename_filt, 'hp')
    preproc.connect(inputnode, ('lowpass', filter_used), rename_filt, 'lp')
    
    """
    Intensity Normalization
    """
    
    # Scale mean value of run to 10000
    meanscale = pe.MapNode(interface=fsl.ImageMaths(suffix='_gms'),
                          iterfield=['in_file', 'op_string'],
                          name='08_meanscale')
    preproc.connect(selectnode, 'out', meanscale, 'in_file')
    ## function to get scaling factor
    preproc.connect(medianval, ('out_stat', getmeanscale), meanscale, 'op_string')
    
    """
    Get Final Functional Data and Mean
    """
    ## functional 4D
    renamer.connect(meanscale, 'out_file', 'func_preproc')
    ## mean
    meanfunc4 = pe.MapNode(fsl.MeanImage(),
                           iterfield=["in_file"],
                           name="08_meanfunc4")
    preproc.connect(meanscale, 'out_file', meanfunc4, 'in_file')
    renamer.connect(meanfunc4, 'out_file', 'func_mean')
    
    """
    Combine all ? (maybe make another workflow)
    """
    
    # Use the -t Merge option to 
    
    # Use fslmaths -Tmin -bin to get proper mask
    
    # Get mean image across entire run
    
    # Combine all the motion together
    
    return preproc


class FuncPreprocParser(usage.NiParser):
    def _create_parser(self, *args, **kwrds):
        """Create command-line interface"""
        parser = super(FuncPreprocParser, self)._create_parser(
            description="""
                Run preprocessing for each participant's functional images.
            """
        )
        group = parser.add_argument_group('Functional Preprocessing Options')
        group.add_argument('--label', required=True)
        group.add_argument('--func', nargs=2, action=usage.store_io, required=True)
        group.add_argument('--fwhm', type=float, default=5.0)
        group.add_argument('--hpfilter', type=float, default=128)
        group.add_argument('--lpfilter', type=float, default=-1)
        group.add_argument('--tr', type=float, required=True)
        group.add_argument('--motion-nstages', type=int, choices=range(1,5), default=3)
        group.add_argument('--whichvol', default="middle", choices=["first", "middle", "mean"])
        group.add_argument("--timeshift", default=False, action="store_true")
        group.add_argument("--tpattern", choices=["alt+z", "alt+z2", "alt-z", "seq+z", "seq-z"])
        return parser
    
    def _post_run(self):
        """Fix permissions of output"""
        outputs = self.args.outputs
        subject_list = self.args.subject_list
        outputs = [ op.join(outputs.basedir, s, outputs.func) for s in subject_list ]
        for output in outputs:
            p = Process("chmod -R 775 %s" % output, to_print=True)
            if p.retcode != 0:
                print 'Error: chmod -R 775 %s' % output
        return
    
    

def main(arglist):
    pp = FuncPreprocParser()
    pp(functional_preprocessing, arglist)

def test_wf():
    arglist = "-s tb3417 -b /Users/zarrar/Projects/tnetworks /Users/zarrar/Projects/tnetworks/output --workingdir /Users/zarrar/Projects/tnetworks/tmp --tr 2 --fwhm 5 --func *_wm_run*.nii.gz wm --plugin Linear"
    main(arglist.split())

if __name__ == "__main__":
    main(sys.argv[1:])

