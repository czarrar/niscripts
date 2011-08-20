import os, yaml
import os.path as op
from nipype.interfaces.base import (traits, TraitedSpec, CommandLineInputSpec, CommandLine, isdefined, File)

from nipype.utils.filemanip import split_filename
from analysis.base import SubjectBase


class FeatSubjectInputSpec(CommandLineInputSpec):
    config_file = File(desc = 'input file to feat_subject.py, yaml configuration file',
                       argstr = "-c %s",
                       mandatory = True, 
                       exists = True)
    subject = traits.Str(desc = "ID for one subject", 
                            argstr = "-s %s", 
                            mandatory = True)
    combine = traits.Bool(desc = "Combine the inputs before creating the fsf", 
                          argstr = "--combine")
    decon = traits.Bool(desc = "Run AFNI's 3dDeconvolve beforehand", 
                          argstr = "--res-decon")
    fsf = traits.Bool(desc = "Create the fsf design file before running the analysis", 
                      argstr = "--fsf")
    feat = traits.Bool(desc = "Run the analysis", 
                       argstr = "--feat")
    regress = traits.Bool(desc = "Regress out stuff using fsl_regfilt", 
                          argstr = "--regress")
    verbose = traits.Bool(desc = "More output", 
                       argstr = "--verbose")
    debug = traits.Bool(desc = "Lots of output", 
                       argstr = "--debug")
    dry_run = traits.Bool(desc = "Don't execute anything just a dry-run.", 
                       argstr = "--dry-run")
    log_dir = traits.Directory(desc = "Log directory", 
                       argstr = "--log-dir")


class FeatSubjectOutputSpec(TraitedSpec):
    func = File(desc="Output from combining multiple functional runs")
    decon_res = File(desc="Output from 3dDeconvolve")
    decon_mat = File(desc="Output from 3dDeconvolve")
    confound = File(desc="Confound matrix to regress out")
    fsf = File(desc="Design file for FEAT analysis")
    mat = File(desc="Mat file generated from fsf file")
    feat = traits.Directory(desc="Output directory of feat analysis")
    regfilt = File(desc="Residuals (with mean) from --regress")

class FeatSubject(CommandLine):
    input_spec = FeatSubjectInputSpec
    output_spec = FeatSubjectOutputSpec
    cmd = 'feat_subject_worker.py'
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        f = file(self.inputs.config, 'r')
        config_dict = yaml.load(f)
        f.close()
        
        template_vars = config_dict.pop("vars", {})
        template_vars['subject'] = self.inputs.subject
        sub = SubjectBase(verbosity=0, template_vars=template_vars)
        
        if 'combine' in config_dict:
            outputs['func'] = sub._substitute(config_dict['combine']['data']['outfunc'])
        
        if 'res_decon' in config_dict:
            if 'outfunc' in config_dict['res_decon']['data']:
                outputs['decon_res'] = sub._substitute(
                                            config_dict['res_decon']['data']['outfunc'])
            if 'outmat' in config_dict['res_decon']['data']:
                outputs['decon_mat'] = sub._substitute(
                                            config_dict['res_decon']['data']['outmat'])
        
        if 'fsf' in config_dict:
            outputs['fsf'] = sub._substitute(config_dict['combine']['data']['outfunc'])
            outputs['mat'] = op.splitext(outputs['fsf'])[0] + ".mat"
            outputs['feat'] = sub._substitute(config_dict['combine']['data']['outdir'])
        
        if 'regress' in config_dict:
            outputs['regfilt'] = sub._substitute(config_dict['combine']['data']['out_file'])
        
        return outputs
    

