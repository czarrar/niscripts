import argparse, os, re
from nipype.interfaces.base import Bunch

#####
# Additional Actions
#####

class store_input(argparse.Action):
    def __init__(self, check_file=False, check_dir=False, **kwrds):
        super(store_input, self).__init__(**kwrds)
        self.check_file = check_file
        self.check_dir = check_dir
    
    def __call__(self, parser, namespace, value, option_string=None):
        if self.check_file and not os.path.isfile(value):
            parser.error("File '%s' from option '%s' does not exist" % (value, option_string))
        if self.check_dir and not os.path.isdir(value):
            parser.error("Dir '%s' from option '%s' does not exist" % (value, option_string))
        namespace.inputs[self.dest] = value
    

class store_output(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        namespace.outputs[self.dest] = value
    

class store_io(argparse.Action):
    def __init__(self, check_file=False, check_dir=False, **kwrds):
        super(store_io, self).__init__(**kwrds)
        self.check_file = check_file
        self.check_dir = check_dir
    
    def __call__(self, parser, namespace, values, option_string=None):
        if self.check_file and not os.path.isfile(values[0]):
            parser.error("File '%s' from option '%s' does not exist" % (values[0], option_string))
        if self.check_dir and not os.path.isdir(values[0]):
            parser.error("Dir '%s' from option '%s' does not exist" % (values[0], option_string))
        namespace.inputs[self.dest] = values[0]
        namespace.outputs[self.dest] = values[1]
    

class store_filename(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if not os.path.isfile(value):
            parser.error("File '%s' does not exist" % value)
        setattr(namespace, self.dest, os.path.abspath(value))
    

class store_directory(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if not os.path.isdir(value):
            parser.error("Directory '%s' does not exist" % value)
        setattr(namespace, self.dest, os.path.abspath(value))
    

class store_plugin(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values[0] == "MultiProc":
            nvalues = 2
            namespace.plugin_args = {'n_procs': int(values[1])}
        else:
            nvalues = 1
            namespace.plugin_args = {}
        if len(values) > nvalues:
            parser.error("Too many arguments (%i) specified for %s" % (len(values), option_string))
        setattr(namespace, self.dest, values[0])
    

#####
# ArgumentParser
#####

class NiArgumentParser(argparse.ArgumentParser):
    __line2args = re.compile("(?P<opt>[\w\-]+)[:]\ *(?P<arg>.*)")
    
    def parse_args(self, args, namespace=None, **kwrds):
        if namespace is None:
            namespace = argparse.Namespace()
        namespace.inputs = {}
        namespace.outputs = {}
        return super(NiArgumentParser, self).parse_args(args, namespace, **kwrds)
    
    def convert_arg_line_to_args(self, arg_line):
        results = self.__line2args.search(arg_line).groupdict()
        if results["opt"]:
            if len(results["opt"]) == 1:
                arg_line = "-%s %s" % (results["opt"], results["arg"])
            elif len(results["opt"]) > 1:
                arg_line = "--%s %s" % (results["opt"], results["arg"])
            else:
                raise Exception("opts is weird: %s" % results["opt"])
        else:
            arg_line = results["arg"]
        for arg in arg_line.split():
            if not arg.strip():
                continue
            yield arg
    

class NiParser(object):
    def __init__(self, *args, **kwrds):
        self.parser = self._create_parser(*args, **kwrds)
        self.workflow = None
    
    def __call__(self, *args, **kwrds):
        return self.run(*args, **kwrds)
    
    def _create_parser(self, fileprefix='@', defaults=argparse.SUPPRESS, **kwrds):
        parser = NiArgumentParser(fromfile_prefix_chars=fileprefix, 
                                  argument_default=defaults, **kwrds)
        group = parser.add_argument_group('Standard Options')
        group.add_argument('-s', '--subjects', nargs="+", dest="subject_list", required=True)
        group.add_argument('-b', '--basedir', nargs=2, action=store_io, required=True, check_dir=True)
        group.add_argument('--workingdir', action=store_directory, required=True)
        group.add_argument('--output-type', default="NIFTI_GZ", choices=['NIFTI', 'NIFTI_GZ'])
        group.add_argument('--name')
        group.add_argument('--plugin', nargs="+", action=store_plugin, required=True)
        return parser
    
    def run(self, procfun, arglist):
        args = self.parser.parse_args(arglist)
        kwrds = vars(args)
        plugin = kwrds.pop('plugin')
        plugin_args = kwrds.pop('plugin_args')
        kwrds["inputs"] = Bunch(**kwrds["inputs"])
        kwrds["outputs"] = Bunch(**kwrds["outputs"])
        workflow = procfun(**kwrds)
        self.workflow = workflow
        return workflow.run(plugin=plugin, plugin_args=plugin_args)
        
