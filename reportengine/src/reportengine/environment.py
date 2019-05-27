# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 00:09:52 2016

@author: Zahari Kassabov
"""
import pathlib
import logging
from collections.abc import Sequence
import shutil

log = logging.getLogger(__name__)

class EnvironmentError_(Exception): pass

available_figure_formats = {
 'eps': 'Encapsulated Postscript',
 'jpeg': 'Joint Photographic Experts Group',
 'jpg': 'Joint Photographic Experts Group',
 'pdf': 'Portable Document Format',
 'pgf': 'PGF code for LaTeX',
 'png': 'Portable Network Graphics',
 'ps': 'Postscript',
 'raw': 'Raw RGBA bitmap',
 'rgba': 'Raw RGBA bitmap',
 'svg': 'Scalable Vector Graphics',
 'svgz': 'Scalable Vector Graphics',
 'tif': 'Tagged Image File Format',
 'tiff': 'Tagged Image File Format'
}

class Environment:
    def __init__(self, *, output=None, formats=('pdf',),
                 default_figure_format=None, loglevel=logging.DEBUG,
                 config_yml = None,
                 **kwargs):
        if output:
            self.output_path = pathlib.Path(output).absolute()
        else:
            self.output_path = output
        self.figure_formats = formats
        self._default_figure_format = default_figure_format
        self.loglevel = loglevel
        self.extra_args = kwargs
        self.config_yml = config_yml

    @property
    def figure_formats(self):
        return self._figure_formats

    @property
    def default_figure_format(self):
        if self._default_figure_format is None:
            return self.figure_formats[0]
        else:
            return self._default_figure_format

    @property
    def config_rel_path(self):
        """A relative path with respect to the config file, or the current
        PWD as a fallback."""
        if self.config_yml:
            return pathlib.Path(self.config_yml).parent
        return pathlib.Path('.')


    @default_figure_format.setter
    def default_figure_format(self, fmt):
        self._default_figure_format = fmt

    @figure_formats.setter
    def figure_formats(self, figure_formats):
        if isinstance(figure_formats, str):
            figure_formats = (figure_formats,)
        if not isinstance(figure_formats, Sequence):
            raise EnvironmentError_("Bad figure format specification: %s. "
                                    "Must be a string or a list." % figure_formats)

        bad_formats = set(figure_formats) - set(available_figure_formats)
        if bad_formats:
            raise EnvironmentError_("The following are not valid figure"
            "formats %s:\nIt must be one of:\n%s" % (bad_formats,
            '\n'.join('%s: %s'%(k,v) for k,v in available_figure_formats.items())))
        self._figure_formats = figure_formats

    def init_output(self):
        if self.output_path and self.output_path.is_dir():
            log.warning("Output folder exists: %s Overwriting contents" %
                     self.output_path)
        else:
            try:
                self.output_path.mkdir()
            except OSError as e:
                raise EnvironmentError_(e) from e
        self.input_folder = self.output_path/'input'
        self.input_folder.mkdir(exist_ok=True)
        if self.config_yml:
            try:
                shutil.copy2(self.config_yml, self.input_folder/'runcard.yaml')
            except shutil.SameFileError:
                pass

        #TODO: Decide if we want to create these always or not
        self.figure_folder = (self.output_path/'figures')
        self.figure_folder.mkdir(exist_ok=True)

        self.table_folder = (self.output_path/'tables')
        self.table_folder.mkdir(exist_ok=True)

    def get_figure_paths(self, handle):
        for fmt in self.figure_formats:
            yield self.figure_folder / (handle + '.' + fmt)

    @classmethod
    def ns_dump_description(cls):
        return dict(
            output_path = "Folder where the results are to be written.",
            config_rel_path = cls.config_rel_path.__doc__,

        )

    def ns_dump(self):
        return {k: getattr(self, k) for k in self.ns_dump_description()
                if getattr(self, k) is not None}
