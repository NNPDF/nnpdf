# -*- coding: utf-8 -*-
from __future__ import print_function
import sys
from setuptools import setup

if sys.version_info < (3,5):
    print("postfit requires Python 3.5 or later", file=sys.stderr)
    sys.exit(1)

with open('README.md') as f:
    long_desc = f.read()

setup(name= "make-fit",
      version = '0.1',
      description = "Filters bad replicas and produces final grid.",
      author = "Zahari Kassabov",
      author_email = "kassabov@to.infn.it",
      url="https://gitlab.cern.ch/NNPDF/nnpdfcpp",
      long_description = long_desc,
      packages = ['makefit'],
      entry_points = {'console_scripts':
                    ['make-fit = makefit:main']},
      zip_safe = False,
      classifiers=[
            'Operating System :: Unix',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics',
            ],
     )

