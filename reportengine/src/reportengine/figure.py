 # -*- coding: utf-8 -*-
"""
Save generated figures in the correct path. Use::

    @figure
    def provider(arg):
       return plt.figure(...)

to have the figure be automatically saved in the correct path, once it is
constructed. Similarly use::

    @figuregen
    def provider(arg):
       for ...:
           yield plt.figure(...)

to have the action applied to each element of a generator.

The figures will be automatically closed.

Created on Thu Mar 10 00:59:31 2016

@author: Zahari Kassabov
"""
import logging

import numpy as np

from reportengine.formattingtools import  spec_to_nice_name
from reportengine.utils import add_highlight, normalize_name

__all__ = ['figure', 'figuregen']

log = logging.getLogger(__name__)

class Figure():
    def __init__(self, paths):
        self.paths = paths

    @property
    def as_markdown(self):
        links = ' '.join('[{ext}]({path})'.format(ext=path.suffix, path=path) for path in self.paths)
        return '![{1}]({0})\n'.format(self.paths[0], links)


def prepare_paths(*,spec, namespace, environment ,**kwargs):
    paths = environment.get_figure_paths(spec_to_nice_name(namespace, spec))
    #list is important here. The generator gives a hard to trace bug when
    #running in parallel
    return {'paths':list(paths), 'output':environment.output_path}

def savefig(fig, *, paths, output ,suffix=''):
    """Final action to save figures, with a nice filename"""
    #Import here to avoid problems with use()
    import matplotlib.pyplot as plt

    outpaths = []
    for path in paths:
        if suffix:
            suffix = normalize_name(suffix)
            path = path.with_name('_'.join((path.stem, suffix)) + path.suffix)
        log.debug("Writing figure file %s" % path)

        #Numpy can produce a lot of warnings while working on producing figures
        with np.errstate(invalid='ignore'):
            fig.tight_layout()
            fig.savefig(str(path))
        outpaths.append(path.relative_to(output))
    plt.close(fig)
    return Figure(outpaths)

def savefiglist(figures, paths, output):
    """Final action to save lists of figures. It adds a numerical index as
    a suffix, for each figure in the generator."""

    res = []

    for i, fig in enumerate(figures):
        #Support tuples with (suffix, figure)
        if isinstance(fig, tuple):
            fig, suffix = fig
        else:
            suffix = str(i)
        res.append(savefig(fig, paths=paths, output=output, suffix=suffix))
    return res


@add_highlight
def figure(f):
    f.prepare = prepare_paths
    f.final_action = savefig
    return f

@add_highlight
def figuregen(f):
    f.prepare = prepare_paths
    f.final_action = savefiglist
    return f
