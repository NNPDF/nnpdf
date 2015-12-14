#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Stefano Carrazza'
__license__ = 'GPL'
__version__ = '1.0.0'
__email__ = 'stefano.carrazza@cern.ch'

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import subprocess
import sys

def call_command(command):
    l = command.split()
    try:
        result = subprocess.check_output(l)
    except OSError as e:
        print("Could not call %s: %s.\n"
              "Please make sure the relevant command is installed."
              % (l[0], e.strerror) )
        sys.exit(1)
    return result.decode().rstrip()

ext = Extension("nnpdf", ['*.pyx'],
                include_dirs = [call_command('lhapdf-config --incdir'),
                                '../src/NNPDF'],
                libraries = ['nnpdf','LHAPDF','gsl'],
                library_dirs = [call_command('nnpdf-config --libdir'),
                                call_command('lhapdf-config --libdir')],
                language = 'c++',
                )

setup(
    name = 'nnpdf',
    author = 'Stefano Carrazza et al.',
    author_email = 'stefano.carrazza@cern.ch',
    version = call_command('nnpdf-config --version'),
    ext_modules = cythonize(ext),
    )
