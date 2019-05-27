# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:42:17 2016

@author: Zahari Kassabov
"""
import traceback
import logging
import copy

import pygments
from pygments import lexers
from pygments.formatters import TerminalFormatter
import blessings



def color_exception(etype, evalue, tb):
    txt = ''.join(traceback.format_exception(etype, evalue, tb))
    lex = lexers.get_lexer_by_name('py3tb')
    tb = pygments.highlight(txt, lex, TerminalFormatter())
    return tb

t = blessings.Terminal()

class ColorHandler(logging.StreamHandler):
    colors = {
        logging.DEBUG: {'[%(levelname)s]:': t.bold},
        logging.INFO: {'[%(levelname)s]:': t.bold_green},
        logging.WARNING: {'[%(levelname)s]:': t.bold_yellow},
        logging.ERROR: {'[%(levelname)s]:': t.bold_red,
                          '%(message)s': t.bold},
        logging.CRITICAL: {'[%(levelname)s]:': t.bold_white_on_red,
                           '%(message)s': t.bold},
    }

    _fmt = '[%(levelname)s]: %(message)s'

    def new_formatter(self, fmt):
        if self.formatter:
            cls = type(self.formatter)
            return cls(fmt, self.formatter.datefmt, self.formatter._style)
        else:
            return logging.Formatter(fmt)

    def color_record_copy(self, record):
        record = copy.copy(record)
        return record

    def setFormatter(self, formatter):
        # HACK: peeping format string passed by user to `logging.Formatter()`
        if formatter._fmt:
            self._fmt = formatter._fmt
        super().setFormatter(formatter)


    def format(self, record):
        levelcolors = self.colors[record.levelno]
        fmt = self._fmt
        for s, subs in levelcolors.items():
            fmt = fmt.replace(s, subs(s))
        return self.new_formatter(fmt).format(record)


