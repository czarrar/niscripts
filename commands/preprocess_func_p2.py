#!/usr/bin/env python
"""
This script uses nipype for functional preprocessing
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

def get_overlay_args(fname):
    """Args for the overlay1 option of slicer"""
    return (fname, 1, 1)

tolist = lambda x: [x]

def functional_preprocessing(
    subject_list, 
    inputs, outputs, workingdir, output_type, 
    fwhm, hpfilter, lpfilter, tr, label, outlabel, 
    name="%s_preprocessing_p2"):
    """todo"""
    
    #####
    # Setup pipeline
    #####
    
    if name.find("%s") != -1:
        name = name % label
    else:
        print "ERROR: You must have a '%s' in the name for the label"
        raise SystemExit(2)
    
    
    func_preproc_name = "func_preproc_%s_%ism_%ihp_%ilp" % (outlabel, int(fwhm*10), 
                            hpfilter*(hpfilter>0), lpfilter*(lpfilter>0))
    
    preproc = create_func_preproc_workflow(name, func_preproc_name)
    preproc.base_dir = workingdir
    
    preproc.inputs.inputspec.fwhm = fwhm    
    preproc.inputs.inputspec.highpass = float(hpfilter)/(2*tr)
    preproc.inputs.inputspec.lowpass = float(lpfilter)/(2*tr)
    
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
                                                    outfields=['func', 'func_mask'],
                                                    sort_filelist=True), 
                         name='datasource')
    datasource.inputs.base_directory=os.path.abspath(inputs.basedir)
    datasource.inputs.template = "*"
    datasource.inputs.field_template = dict(
        func = os.path.join("%s", inputs.funcdir, 'run_?', inputs.infunc),
        func_mask = os.path.join("%s", inputs.funcdir, 'run_?', inputs.inmask)
    )
    datasource.inputs.template_args = dict(
                                        func = [['subject_id']],
                                        func_mask = [['subject_id']]
                                      )
    
    # Get the number of runs for each participant (for renaming purposes)
    datasource.inputs.subject_id = subject_list
    ds = datasource.run()
    print [ x for x in ds.outputs.func ]
    print ds.outputs.func
    raise SystemExit(1)
    runs = [ len(x) for x in ds.outputs.func ]  # list with number of runs per subject
    max_runs = max(runs)
    
    # Link inputs
    preproc.connect(subinfo, 'subject_id', datasource, 'subject_id')
    preproc.connect(datasource, 'func', inputnode, 'func')
    preproc.connect(datasource, 'func_mask', inputnode, 'func_mask')
    
    
    ######
    # Setup data sink
    ######
    
    # Datasink
    ## will get: "base_directory/subject_id/output_anatdir"
    datasink = pe.Node(interface=nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = outputs.basedir
    
    ## substitute map stuff
    preprocsinksubs = get_mapnode_substitutions(preproc, outputnode, max_runs)
    datasink.inputs.substitutions = preprocsinksubs
    
    # replace subject_id stuff with functional scan
    datasink.inputs.regexp_substitutions = (r"_subject_id_(\w|\d)+", outputs.funcdir)
    
    ## connect
    preproc.connect(subinfo, 'subject_id', datasink, 'container')
    output_fields = outputnode.outputs.get()
    for field in output_fields:
        preproc.connect(outputnode, field, datasink, "@%s" % field)
    
    return preproc

# whichvol can be first, middle, or mean
# hp should be in TRs
def create_func_preproc_workflow(name='functional_preprocessing', func_preproc_name="func_preproc"):
        
    #####
    # Setup workflow
    #####
    
    preproc = pe.Workflow(name=name)
    
    
    #####
    # Setup input node
    #####
    
    input_fields = [
        "func",         # list of functional images to run
        "func_mask",    # list of functional images to run
        "fwhm",         # fwhm
        "highpass",     # highpass filtering
        "lowpass",      # lowpass filtering
    ]
    inputnode = pe.Node(interface=util.IdentityInterface(fields=input_fields), 
                            name="inputspec")
    
    
    #####
    # Setup output node
    #####
    
    # Outputs
    output_fields = [
        func_preproc_name,  # final output (also has been intensity normalized)
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
    
    """
    Smoothing
    """
    
    # Anisitropic-like? smoothing with your friend susan
    smooth = create_susan_smooth(name="01_susan_smooth")
    preproc.connect(inputnode, 'fwhm', smooth, 'inputnode.fwhm')
    preproc.connect(img2float, 'out_file', smooth, 'inputnode.in_files')
    preproc.connect(inputnode, 'func_mask', smooth, 'inputnode.mask_file')
    
    # Mask smoothed data with dilated mask
    funcbrain = pe.MapNode(fsl.ApplyMask(),
                           iterfield=["in_file", "mask_file"],
                           name="02_funcbrain")
    preproc.connect(smooth, 'outputnode.smoothed_files', funcbrain, 'in_file')
    preproc.connect(inputnode, 'func_mask', funcbrain, 'mask_file')
    
    # Determine if want to take output from smoothed or non-smoothed data forward
    # This would be the case if the fwhm is less than 1/3 of the voxel size
    ## gather 2 types of functional images
    concatnode = pe.Node(interface=util.Merge(2),
                         name='02_concat')
    preproc.connect(img2float, ('out_file', tolist), concatnode, 'in1')
    preproc.connect(funcbrain, ('out_file', tolist), concatnode, 'in2')
    ## select one
    selectnode = pe.Node(interface=util.Select(),name='02_select')
    preproc.connect(concatnode, 'out', selectnode, 'inlist')
    preproc.connect(inputnode, ('fwhm', chooseindex), selectnode, 'index')
    
    """
    Filter
    """
    
    filt = pe.MapNode(fsl.TemporalFilter(),
                          iterfield=["in_file"],
                          name="03_filter")
    preproc.connect(inputnode, 'highpass', filt, 'highpass_sigma')
    preproc.connect(inputnode, 'lowpass', filt, 'lowpass_sigma')
    preproc.connect(selectnode, 'out', filt, 'in_file')
    
    """
    Intensity Normalization
    """
    
    # Determine the median value of the functional runs using the mask
    medianval = pe.MapNode(fsl.ImageStats(op_string="-k %s -p 50"),
                           iterfield = ["in_file", "mask_file"],
                           name="04_medianval")
    preproc.connect(filt, 'out_file', medianval, 'in_file')
    preproc.connect(inputnode, 'func_mask', medianval, 'mask_file')
    
    # Scale mean value of run to 10000
    meanscale = pe.MapNode(interface=fsl.ImageMaths(suffix='_gms'),
                          iterfield=['in_file', 'op_string'],
                          name='04_meanscale')
    preproc.connect(filt, 'out_file', meanscale, 'in_file')
    ## function to get scaling factor
    preproc.connect(medianval, ('out_stat', getmeanscale), meanscale, 'op_string')
    
    """
    Get Final Functional Data and Mean
    """
    ## functional 4D
    renamer.connect(meanscale, 'out_file', func_preproc_name)
    
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
        group.add_argument('--outlabel', required=True)
        group.add_argument('--funcdir', nargs=2, action=usage.store_io, required=True)
        group.add_argument('--infunc', action=usage.store_input, required=True)
        group.add_argument('--inmask', action=usage.store_input, required=True)
        group.add_argument('--fwhm', type=float, default=5.0)
        group.add_argument('--hpfilter', type=float, default=128)
        group.add_argument('--lpfilter', type=float, default=-1)
        group.add_argument('--tr', type=float, required=True)
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

