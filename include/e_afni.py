from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.afni.base import AFNITraitedSpec, AFNICommand, Info
from nipype.interfaces.base import Bunch, TraitedSpec, File, Directory, traits, isdefined

import os

class AFNICommandGenFile(AFNICommand):
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
            Flag to change the filename extension to the AFNI output type.
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
        
    def _gen_filename(self, name):
        if name == 'out_file':
            return self._gen_outfilename()
        return None
    

class ThreedcopyInputSpec(AFNITraitedSpec):
    in_file = File(desc = 'input file to 3dcopy',
                    argstr = '%s',
                    position = -2,
                    mandatory = True,
                    exists = True)
    out_file = File(desc = 'output file from 3dcopy',
                    argstr = '%s',
                    position = -1,
                    genfile = True)

class ThreedcopyOutputSpec(AFNITraitedSpec):
    out_file = File(desc = 'copied file',
                    exists = True)

class Threedcopy(AFNICommandGenFile):
    """Copies an image of one type to an image of the same or different type
    using 3dcopy command."""
    
    _cmd = '3dcopy'
    input_spec = ThreedcopyInputSpec
    output_spec = ThreedcopyOutputSpec
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = self._gen_outfilename()
        return outputs
    
    def _gen_outfilename(self):
        out_file = self.inputs.out_file
        if not isdefined(out_file) and isdefined(self.inputs.in_file):
            out_file = self._gen_fname(self.inputs.in_file, suffix="_copy")
        return out_file   
    

class ThreedresampleInputSpec(AFNITraitedSpec):
    in_file = File(desc = 'input file to 3dresample',
                  argstr = '-inset %s',
                  position = -1,
                  mandatory = True,
                  exists = True)
    out_file = File(desc = 'output file from 3dresample',
                   argstr = '-prefix %s',
                   position = -2,
                   genfile = True, xor=['out_prefix'])
    orientation = traits.Str(desc = 'new orientation code',
                             argstr = '-orient %s')
    out_prefix = traits.Str(desc = 'filename of output without path or extension', 
                            xor=['out_file'])
    
class ThreedresampleOutputSpec(AFNITraitedSpec):
    out_file = File(desc = 'reoriented or resampled file',
                    exists = True)

class Threedresample(AFNICommandGenFile):
    """Resample or reorient an image using AFNI 3dresample command.
    
    For complete details, see the `3dresample Documentation.
    <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dresample.html>`_
    """
    
    _cmd = '3dresample'
    input_spec = ThreedresampleInputSpec
    output_spec = ThreedresampleOutputSpec
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = self._gen_outfilename()
        return outputs
    
    def _gen_outfilename(self):
        out_file = self.inputs.out_file
        if not isdefined(out_file) and isdefined(self.inputs.in_file):
            if isdefined(self.inputs.out_prefix):
                pre_out_file = os.path.join(
                                    os.path.dirname(self.inputs.in_file), 
                                    self.inputs.out_prefix
                                )
                suffix = ''
            else:
                pre_out_file = self.inputs.in_file
                suffix = []
                if self.inputs.orientation:
                    suffix.append("_reorient")
                suffix = "+".join(suffix)
            out_file = self._gen_fname(pre_out_file, suffix=suffix)
        return out_file
    

class ThreedSkullStripInputSpec(AFNITraitedSpec):
    in_file = File(desc = 'input file to 3dSkullStrip',
                    argstr = '-input %s',
                    position = -2,
                    mandatory = True,
                    exists = True)
    out_file = File(desc = 'output file from 3dSkullStrip',
                    argstr = '-prefix %s',
                    position = -1,
                    genfile = True)
    options = traits.Str(desc = 'additional options',
                         argstr = '%s')


class ThreedSkullStripOutputSpec(AFNITraitedSpec):
    out_file = File(desc = 'skull stripped image',
                    exists = True)


class ThreedSkullStrip(AFNICommandGenFile):
    """Copies an image of one type to an image of the same or different type
    using 3dcopy command."""
    
    _cmd = '3dSkullStrip'
    input_spec = ThreedSkullStripInputSpec
    output_spec = ThreedSkullStripOutputSpec
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = self._gen_outfilename()
        return outputs
    
    def _gen_outfilename(self):
        out_file = self.inputs.out_file
        if not isdefined(out_file) and isdefined(self.inputs.in_file):
            out_file = self._gen_fname(self.inputs.in_file, suffix="_brain")
        return out_file
    

