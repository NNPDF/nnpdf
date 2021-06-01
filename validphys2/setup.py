from __future__ import print_function
import sys
from setuptools import setup, find_packages

if sys.version_info < (3,6):
    print("validphys requires Python 3.6 or later", file=sys.stderr)
    sys.exit(1)

with open('README.md') as f:
    long_desc = f.read()

setup(name= "validphys",
      version = '4.0',
      description = "NNPDF analysis framework",
      author = "Zahari Kassabov",
      author_email = "kassabov@to.infn.it",
      url="https://gitlab.cern.ch/NNPDF/validphys2",
      long_description = long_desc,
      entry_points = {'console_scripts':
                    [
                        'validphys = validphys.scripts.main:main',
                        'vp-upload = validphys.scripts.vp_upload:main',
                        'wiki-upload = validphys.scripts.wiki_upload:main',
                        'postfit = validphys.scripts.postfit:main',
                        'vp-get = validphys.scripts.vp_get:main',
                        'vp-comparefits = validphys.scripts.vp_comparefits:main',
                        'vp-fitrename = validphys.scripts.vp_fitrename:main',
                        'vp-checktheory = validphys.scripts.vp_checktheory:main',
                        'vp-rebuild-data = validphys.scripts.vp_rebuild_data:main',
                        'vp-pdfrename = validphys.scripts.vp_pdfrename:main',
                        'vp-pdffromreplicas = validphys.scripts.vp_pdffromreplicas:main',
                        'vp-list = validphys.scripts.vp_list:main',
                        'vp-nextfitruncard = validphys.scripts.vp_nextfitruncard:main',
                        'vp-hyperoptplot = validphys.scripts.vp_hyperoptplot:main',
                        'vp-deltachi2 = validphys.scripts.vp_deltachi2:main',
                    ]},
      package_dir = {'': 'src'},
      packages = find_packages('src'),
      package_data = {
           #TODO: Get rid of this nonsense
            '':['*.template', '*.mplstyle', '*.csv', '*.yaml', '*.md', '*.png', '*.pickle'],
            'tests/regressions': ['*'],
            'comparefit': ['*'],
            'tests/baseline': ['*'],
            'cuts': ['*'],
            'scalevariations': ['*'],
            'hyperplot': ['*'],
            'deltachi2': ['*'],
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

