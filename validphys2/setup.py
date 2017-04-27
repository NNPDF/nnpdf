from __future__ import print_function
import sys
from setuptools import setup, find_packages

if sys.version_info < (3,6):
    print("validphys requires Python 3.6 or later", file=sys.stderr)
    sys.exit(1)

with open('README.md') as f:
    long_desc = f.read()

setup(name= "validphys",
      version = '2.0b2',
      description = "NNPDF analysis framework",
      author = "Zahari Kassabov",
      author_email = "kassabov@to.infn.it",
      url="https://gitlab.cern.ch/NNPDF/validphys2",
      long_description = long_desc,
      entry_points = {'console_scripts':
                    [
                        'validphys = validphys.scripts.main:main',
                        'vp-upload = validphys.scripts.vp_upload:main',
                    ]},
      package_dir = {'': 'src'},
      packages = find_packages('src'),
       package_data = {
            '':['*.template', '*.mplstyle']
       },
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

