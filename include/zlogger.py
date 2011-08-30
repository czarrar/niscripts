#!/usr/bin/env python

import os, logging, sys
import os.path as op
import logging.config
from execute import Process

# Add level names
# warning = 30
logging.IMPORTANT = 28  
logging.TITLE = 26
logging.SUBTITLE = 24
logging.COMMAND = 18
logging.COMMAND_INFO = 15
logging.COMMAND_STDOUT = 16
logging.COMMAND_STDERR = 17
logging.addLevelName(logging.IMPORTANT, 'IMPORTANT')
logging.addLevelName(logging.TITLE, 'TITLE')
logging.addLevelName(logging.SUBTITLE, 'SUBTITLE')
logging.addLevelName(logging.COMMAND, 'COMMAND')
logging.addLevelName(logging.COMMAND_INFO, 'COMMAND_INFO')
logging.addLevelName(logging.COMMAND_STDOUT, 'COMMAND_STDOUT')
logging.addLevelName(logging.COMMAND_STDERR, 'COMMAND_STDERR')

# Get logger wrapper
def getLogger(name=None, level=logging.DEBUG, fname=None, extra={}, allow_exceptions=False):
    # create logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.propagate = 0
    
    # create file handler
    if fname is not None:
        fh = logging.FileHandler(fname)
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    
    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    formatter = ColoredFormatter("%(levelcolor)s%(message)s%(endcolor)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    return Adapter(logger, extra, allow_exceptions)

class LoggerError(Exception):
    """todo"""

class LoggerCritical(Exception):
    """todo"""

class Adapter(logging.LoggerAdapter):
    def __init__(self, logger, extra, allow_exceptions=False, **kwargs):
        super(Adapter, self).__init__(logger, extra, **kwargs)
        self.allow_exceptions = allow_exceptions
    
    def command(self, cmd, cwd=None, shell=False, *args, **kwargs):
        self.log(logging.COMMAND, cmd, *args, **kwargs)
        p = Process(cmd, to_print=False, cwd=cwd, shell=shell)
        if p.stdout:
            self.log(logging.COMMAND_STDOUT, p.stdout)
        if p.stderr:
            self.log(logging.COMMAND_STDERR, p.stderr)
        if self.allow_exceptions and p.retcode != 0:
            raise LoggerError("Error calling: %s" % cmd)
        return p
    
    def drycommand(self, msg, *args, **kwargs):
        self.log(logging.COMMAND, msg, *args, **kwargs)
    
    def subtitle(self, msg, *args, **kwargs):
        self.log(logging.SUBTITLE, msg, *args, **kwargs)
    
    def title(self, msg, *args, **kwargs):
        self.log(logging.TITLE, msg, *args, **kwargs)
    
    def important(self, msg, *args, **kwargs):
        self.log(logging.IMPORTANT, msg, *args, **kwargs)
    
    def error(self, msg, *args, **kwargs):
        self.log(logging.ERROR, msg, *args, **kwargs)
        if self.allow_exceptions:
            raise LoggerError(msg)
    
    def critical(self, msg, *args, **kwargs):
        self.log(logging.CRITICAL, msg, *args, **kwargs)
        if self.allow_exceptions:
            raise LoggerCritical(msg)
    
    def fatal(self, msg, *args, **kwargs):
        self.log(logging.CRITICAL, msg, exc_info=1, *args, **kwargs)
        raise SystemExit(1)
    

# Format output (prettify display)
class ColoredFormatter(logging.Formatter):
    """
    Used with logger, formats logger output to have colors on the terminal.
    """
    
    _ansi_code_conversions = {
        'bold': 1,
        'dark': 2,
        'italics': 3,
        'underline': 4,
        'blink': 5,
        'reverse': 7,
        'concealed': 8,
        'cross-out': 9,
        'black': 30,
        'red': 31,
        'green': 32,
        'yellow': 33,
        'blue': 34,
        'magenta': 35,
        'cyan': 36,
        'white': 37,
        'default': 39,
        # 'on_*' convention is 'borrowed' from python package termcolor
        # 'on' indicates the background color
        'on_black': 40,
        'on_red': 41,
        'on_green': 42,
        'on_yellow': 43,
        'on_blue': 44,
        'on_magenta': 45,
        'on_cyan': 46,
        'on_white': 47,
        'on_default': 49
    }
    # list of format options will be converted to comma-seperated numbers with constructor
    _hlevelcolors = {
        'DEBUG': ['magenta'],
        'COMMAND_STDOUT': ['default'],
        'COMMAND_STDERR': ['default'],
        'COMMAND_INFO': ['green', 'dark', 'underline'],
        'COMMAND': ['green', 'bold'],
        'INFO': ['cyan'],
        'SUBTITLE': ['white', 'on_magenta'],
        'TITLE': ['white', 'bold', 'on_blue'],
        'IMPORTANT': ['reverse'],
        'WARNING': ['red'],
        'ERROR': ['white', 'on_red', 'bold'],
        'CRITICAL': ['red', 'on_yellow', 'bold', 'blink'],
    }
    def __init__(self, *args, **kwargs):
        logging.Formatter.__init__(self, *args, **kwargs)
        self._levelcolors = {}
        for k,v in self._hlevelcolors.iteritems():
            self._levelcolors[k] = ";".join([ str(self._ansi_code_conversions[x]) for x in v ])
    def format(self, record):
        record.levelcolor = '\033[' + self._levelcolors[record.levelname] + 'm'
        record.endcolor = '\033[0m'
        return logging.Formatter.format(self, record)

