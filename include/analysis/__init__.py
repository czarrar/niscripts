import os.path as op
import yaml
from copy import deepcopy
from analysis.feat import (CombineSubject, FsfSubject, FsfInfo, FeatSubject)
from analysis.afni import (DeconSubject, RemlSubject, BetaSeriesSubject)

def fromYamlSubject(inputs, run_keys, verbosity=0, dry_run=False, log=None, **user_template_vars):
    """
    Runs both CombineFuncs and FsfSubject
    """
    config = inputs['config']
    if isinstance(config, str):
        if not op.isfile(config):
            raise Exception("Cannot find config file %s" % config)
        config = file(config, 'r')
    config_dict = yaml.load(config)
    
    template_vars = config_dict.pop("vars", {})
    template_vars.update(user_template_vars)
    
    for k in run_keys:
        opts = config_dict.pop(k, None)
        if opts:
            class_name = "".join([ x.capitalize() for x in k.split("_") ]) + "Subject"
            c = eval(class_name)(verbosity, deepcopy(template_vars), dry_run, log)
            c.fromDict(opts)
            c.run()
    
    return

