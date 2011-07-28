import argparse, os

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


class store_basedir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) > 2:
            parser.error("Too many arguments (%i) specified for %s" % (len(values), option_string))
        for value in values:
            if not os.path.isdir(value):
                parser.error("Base directory '%s' does not exist" % value)
        if len(values) == 1:
            values.append(values[0])
        setattr(namespace, 'input_basedir', os.path.abspath(values[0]))
        setattr(namespace, 'output_basedir', os.path.abspath(values[1]))
    


# Parent Parser

parent_parser = argparse.ArgumentParser(add_help=False)

parent_group = parent_parser.add_argument_group('Standard Options')

parent_group.add_argument('-s', '--subjects', nargs="+", required=True, dest="subject_list")

parent_group.add_argument('-b', '--basedir', nargs="+", action=store_basedir, default=argparse.SUPPRESS, required=True)

parent_group.add_argument('--workingdir', action=store_directory, required=True)

parent_group.add_argument('--output-type', default="NIFTI_GZ", choices=['NIFTI', 'NIFTI_GZ'])

parent_group.add_argument('--name', default=argparse.SUPPRESS)
