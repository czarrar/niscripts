from nipype.interfaces.base import (
    traits,
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    isdefined
)

import os
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.fsl.base import Info # todo: make my own Info thing!

class SlicerInputSpec(CommandLineInputSpec):
    in_file = traits.File(desc = 'input file to slicer.py',
                    argstr = '%s',
                    position = -2,
                    mandatory = True,
                    exists = True)
    out_file = traits.File(desc = 'output file from slicer.py',
                    argstr = '%s',
                    position = -1,
                    genfile = True)
    edge_overlay = traits.File(desc = "an outline of a brain image is overlaid on the input", 
                        argstr = '--edge-overlay %s', exists = True)
    overlay1 = traits.Tuple(traits.File, traits.Float, traits.Float, exists = True, 
                            argstr='--overlay %s %.4f %.4f', desc="""image to overlay on input 
                            according to min and max""")
    overlay2 = traits.Tuple(traits.File, traits.Float, traits.Float, exists=True, 
                            argstr='--overlay2 %s %.4f %.4f', xor=['show_negative'], 
                            desc="2nd image to overlay on input according to min and max")
    show_negative = traits.Bool(desc='display negative statistics in overlay', exists = True, 
                                argstr='--show-negative', xor=['overlay2'])
    transparency = traits.Bool(desc='make overlay colors semi-transparent', argstr='-t')
    checkerboard = traits.Bool(desc="use checkerboard mask for overlay", argstr='--checkerboard')
    bg_range = traits.Tuple(traits.Float, traits.Float, argstr='--background-range %.4f %.4f', 
                            desc="setting range of underlay (default: automatic estimation)")
    slice_labels = traits.Bool(desc="display slice numbers", argstr='--slice-labels')
    colour_map = traits.File(desc="use different colour map from that stored in nifti header", 
                            argstr='--lut %s', exists=True)
    scaling = traits.Int(desc="relative size of each slice", argstr="--scale %d")
    intensity_range = traits.Tuple(traits.Float, traits.Float, argstr="--intensity %.3f %.3f", 
                                    desc="min and max intensities to display")
    hide_lr = traits.Bool(desc="Hide (L)eft/(R)ight hemisphere label", argstr="--no-lr-labels")
    slice_name = traits.Enum('a', 'axial', 'c', 'coronal', 's', 'sagittal', argstr="-s %s", 
                            desc="name of slice (a/axial, c/coronal, or s/sagittal)", 
                            mandatory=True)
    width = traits.Int(desc='image width in # of slices', argstr="-w %d", mandatory=True)
    height = traits.Int(desc="image height in # of slices", argstr="-l %d", xor=['sample'])
    sample = traits.Int(desc="include every X # of slices", argstr="-e %d", xor=['height'])
    force = traits.Bool(desc="overwrite output if it exists", argstr="--force")
    verbose = traits.Bool(desc="verbose", argstr="-v")

class SlicerOutputSpec(TraitedSpec):
    out_file = traits.File(desc="PNG image", exists=True)
    
class Slicer(CommandLine):
    input_spec = SlicerInputSpec
    output_spec = SlicerOutputSpec
    cmd = 'slicer.py'
    
    def _gen_fname(self, basename, cwd=None, suffix=None, change_ext=True, ext=None):
        """Generate a filename based on the given parameters.

        The filename will take the form: cwd/basename<suffix><ext>.
        If change_ext is True, it will use the extensions specified in
        <instance>inputs.outputtype.

        Parameters
        ----------
        basename : str
            Filename to base the new filename on.
        cwd : str
            Path to prefix to the new filename. (default is os.getcwd())
        suffix : str
            Suffix to add to the `basename`.  (default is '')
        change_ext : bool
            Flag to change the filename extension to the FSL output type.
            (default True)

        Returns
        -------
        fname : str
            New filename based on given parameters.

        """

        if basename == '':
            msg = 'Unable to generate filename for command %s. ' % self.cmd
            msg += 'basename is not set!'
            raise ValueError(msg)
        if cwd is None:
            cwd = os.getcwd()
        if ext is None:
            ext = Info.outputtype_to_ext(self.inputs.outputtype)
        if change_ext:
            if suffix:
                suffix = ''.join((suffix, ext))
            else:
                suffix = ext
        if suffix is None:
            suffix = ''
        fname = fname_presuffix(basename, suffix = suffix,
                                use_ext = False, newpath = cwd)
        return fname
    
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = self._gen_outfilename()
        return outputs
    
    def _gen_outfilename(self):
        out_file = self.inputs.out_file
        if not isdefined(out_file):
            if isdefined(self.inputs.overlay1):
                self.inputs.overlay1
                tmp_fname, _, _ = self.inputs.overlay1
                out_file = self._gen_fname(tmp_fname, ext='.png')
            elif isdefined(self.inputs.in_file):
                out_file = self._gen_fname(self.inputs.in_file, ext='.png')
        else:
            outputs['out_file'] = self._gen_fname(self.inputs.out_file, ext='.png')
        return out_file
    
    def _gen_filename(self, name):
        if name == 'out_file':
            return self._gen_outfilename()
        return None
    
