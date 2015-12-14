from __future__ import print_function
import sys
import glob
from setuptools import setup, Extension, find_packages
import subprocess
import platform

if sys.version_info < (3,4):
    print("Wrapper requires Python 3.4 or later", file=sys.stderr)
    sys.exit(1)

def call_command(command):
    l = command.split()
    try:
        result = subprocess.check_output(l)
    except OSError as e:
        print("Could not call %s: %s.\n"
              "Please make sure the relevant command is isntalled."
              % (l[0], e.strerror), file=sys.stderr)
        sys.exit(1)
    return result.decode().rstrip()

nnpdf_includes = call_command('nnpdf-config --cppflags').split()
nnpdf_libs = call_command('nnpdf-config --ldflags').split()

extra_compile_args = nnpdf_includes
extra_link_args = nnpdf_libs

if platform.system() == 'Darwin':
    mac_ver = platform.mac_ver()[0]
    extra_compile_args += ['-mmacosx-version-min=%s' % mac_ver]
    extra_link_args += ['-mmacosx-version-min=%s' % mac_ver]

interfaces = glob.glob1("src/wrapper/", "*.i")
names = [i[:-2] for i in interfaces]
sources = ["src/wrapper/_sources/%s_wrap.cxx"%n for n  in names]

ext_modules = [Extension('NNPDF._' + name,
                    extra_compile_args = extra_compile_args ,
                    extra_link_args = extra_link_args,
                    sources = [source], language="c++") 
                    for name,source in zip(names, sources) ]

setup (name = 'NNPDF',
       version = '0.1',
       description = "PDFs for phenomenology",
       author = 'NNPDF Collaboration',
       author_email = 'kassabov@to.infn.it',
       url = 'https://nnpdf.hepforge.org/',
       long_description = "See README",
       ext_modules = ext_modules,
       #http://stackoverflow.com/questions/24799146/use-multiprocessing-process-inside-a-script-installed-by-setuptools
       package_dir = {'': 'src/wrapper'},
       packages = find_packages('src/wrapper'),
       zip_safe = False,
       classifiers=[
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.4',
            ],
       )
