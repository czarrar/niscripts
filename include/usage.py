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
        elif values[0] == "SGE":
            nvalues = 2
            namespace.plugin_args = {'qsub_args': str(values[1])}
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
    _add_inputs = True
    _add_outputs = True
    
    def parse_args(self, args, namespace=None, **kwrds):
        if namespace is None:
            namespace = argparse.Namespace()
        if self._add_inputs:
            namespace.inputs = {}
        if self._add_outputs:
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
        # below created/set by compile
        self.args = None
        self.kwrds = None
        self.plugin = None
        self.plugin_args = None
        self._is_compiled = False
        # below created by run
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
        group.add_argument('--plugin', nargs="+", action=store_plugin, default="Linear")
        group.add_argument('--log-dir', dest="log_dir")
        group.add_argument('--crash-dir', dest="crash_dir")
        return parser
    
    def _pre_compile(self):
        pass
    
    def _post_compile(self):
        pass
    
    def compile(self, arglist, namespace=None):
        self._pre_compile()
        args = self.parser.parse_args(arglist, namespace=namespace)
        self.args = args
        
        kwrds = vars(args)
        self.plugin = kwrds.pop('plugin')
        self.plugin_args = kwrds.pop('plugin_args', {})
        self.log_dir = kwrds.pop('log_dir', None)
        self.crash_dir = kwrds.pop('crash_dir', None)
        if self.log_dir:
            self.log_dir = os.path.abspath(os.path.expanduser(self.log_dir))
        if 'inputs' in kwrds:
            kwrds["inputs"] = Bunch(**kwrds["inputs"])
        if 'outputs' in kwrds:
            kwrds["outputs"] = Bunch(**kwrds["outputs"])
        self.kwrds = kwrds
        
        self._post_compile()
        
        self._is_compiled = True
        return
    
    def _pre_run(self):
        return
    
    def _post_run(self):
        return
    
    def run(self, procfun, arglist=None):
        if arglist is not None:
            self.compile(arglist)
        if not self._is_compiled:
            raise Exception('Called run for NiParser without calling compile or giving arglist')
        
        self._pre_run()
        wf = procfun(**self.kwrds)
        if not isinstance(wf, list):
            wf = [wf]
        self.workflow = wf
        if self.log_dir:
            if not os.path.isdir(self.log_dir):
                os.mkdir(self.log_dir)
            for wf in self.workflow
                wf['config']['log_directory'] = self.log_dir
        if self.crash_dir:
            crashdir = os.path.abspath(os.path.expanduser(self.crash_dir))
            if not os.path.isdir(crashdir):
                os.makedirs(crashdir)
            for wf in self.workflow:
                wf.config = dict(crashdump_dir=crashdir)
        res = [ wf.run(plugin=self.plugin, plugin_args=self.plugin_args) for wf in self.workflow ]
        self._post_run()
        
        return res
    

