project(
'libnnpdf',
'cpp',
version: '1.1',
license: 'SECRET',
default_options : ['cpp_std=c++14']
)

dependencies = [
    dependency('lhapdf'),
    dependency('gsl'),
    dependency('sqlite3'),
    dependency('libarchive'),
    dependency('yaml-cpp'),
]

headers = [
    'src/NNPDF/commondata.h',
    'src/NNPDF/common.h',
    'src/NNPDF/chisquared.h',
    'src/NNPDF/dataset.h',
    'src/NNPDF/exceptions.h',
    'src/NNPDF/experiments.h',
    'src/NNPDF/fastkernel.h',
    'src/NNPDF/fkgenerator.h',
    'src/NNPDF/fkset.h',
    'src/NNPDF/lhapdfset.h',
    'src/NNPDF/logger.h',
    'src/NNPDF/nnmpi.h',
    'src/NNPDF/nnpdfdb.h',
    'src/NNPDF/parametrisation.h',
    'src/NNPDF/pdfset.h',
    'src/NNPDF/positivity.h',
    'src/NNPDF/randomgenerator.h',
    'src/NNPDF/thpredictions.h',
    'src/NNPDF/timer.h',
    'src/NNPDF/utils.h',
    'src/NNPDF/pathlib.h',
]

#handle configure file
conf = configuration_data()
conf.set_quoted('VERSION', meson.project_version())
configure_file(
    input: 'src/NNPDF/config.h.in',
    output: 'config.h',
    configuration: conf
)

#handle pathlib

profile_name = 'nnprofile.yaml'
profile_input = '.'.join([profile_name, 'in'])
profile_conf = configuration_data()

profile_path = join_paths([
   get_option('prefix'),
   'share', 'NNPDF',
])

if get_option('bureocratic_paths')
	profile_conf.set('PROFILE_PREFIX', profile_path)
else
	profile_conf.set('PROFILE_PREFIX', join_paths(['..', 'nnpdfcpp']))
endif

profile_source = configure_file(
    input: profile_input,
	output: profile_name,
    configuration: profile_conf
)

install_data(profile_source, install_dir:profile_path)
add_global_arguments(
    '-DDEFAULT_NNPDF_PROFILE_PATH="@0@"'.format(
	    join_paths([profile_path,profile_name])
    ),
    language:'cpp'
)

sources = [
    'src/common.cc',
    'src/commondata.cc',
    'src/chisquared.cc',
    'src/dataset.cc',
    'src/experiments.cc',
    'src/fastkernel.cc',
    'src/fkconvolute.cc',
    'src/fkgenerator.cc',
    'src/fkset.cc',
    'src/lhapdfset.cc',
    'src/logger.cc',
    'src/nnmpi.cc',
    'src/nnpdfdb.cc',
    'src/parametrisation.cc',
    'src/pdfset.cc',
    'src/positivity.cc',
    'src/randomgenerator.cc',
    'src/thpredictions.cc',
    'src/utils.cc',
    'src/pathlib.cc',
]

install_headers(headers, subdir:'NNPDF')

inc = include_directories(['src', 'src/NNPDF'])

lib = shared_library(
    'nnpdf',
	sources,
	include_directories: inc,
	dependencies: dependencies,
	install: true,
    install_dir: 'lib/',
)

test_sources = [
    'tests/main.cc',
	'tests/experiment_test.cc',
]

run_test = executable(
    'run_test',
	test_sources,
	link_with:[lib],
    include_directories: inc,
	build_by_default: false,
)

test('run_test', run_test)



pkg_mod = import('pkgconfig')
pkg_mod.generate(
    libraries: lib,
	version: meson.project_version(),
	name: 'libnnpdf',
	filebase: 'nnpdf',
	description: 'NNPDF C++ base library.',

)

swig = find_program('swig')
swig_source = 'wrapper/src/nnpdf.i'

swig_gen = generator(swig,
    output: ['@BASENAME@_wrap.cxx'],
	arguments: ['-c++', '-python', '-py3', '-o', './@OUTPUT@',
	    '-outdir', '.', '-I../src', '-MMD', '@INPUT@'],
    depfile: '@BASENAME@_wrap.d',
)

wrapper_cxx = swig_gen.process(swig_source)


wrapper_deps = [
     dependency('lhapdf')
]

python = import('python3')

py3_dep = dependency('python3')

python_install_dir = join_paths(['/', python.sysconfig_path('platlib'), 'NNPDF'])

#TODO: Find a batter way to do this
#Except for __init__.py, these are all deprecated anyway
pyfiles = [
    'wrapper/src/__init__.py',
    'wrapper/src/common.py',
    'wrapper/src/commondata.py',
    'wrapper/src/dataset.py',
    'wrapper/src/experiments.py',
    'wrapper/src/fastkernel.py',
    'wrapper/src/fkgenerator.py',
    'wrapper/src/fkset.py',
    'wrapper/src/lhapdfset.py',
    'wrapper/src/logger.py',
    'wrapper/src/nnmpi.py',
    'wrapper/src/nnpdfdb.py',
    'wrapper/src/parametrisation.py',
    'wrapper/src/pathlib.py',
    'wrapper/src/pdfset.py',
    'wrapper/src/positivity.py',
    'wrapper/src/randomgenerator.py',
    'wrapper/src/thpredictions.py',
    'wrapper/src/timer.py',
    'wrapper/src/utils.py',


    'wrapper/src/_common.py',
    'wrapper/src/_commondata.py',
    'wrapper/src/_dataset.py',
    'wrapper/src/_experiments.py',
    'wrapper/src/_fastkernel.py',
    'wrapper/src/_fkgenerator.py',
    'wrapper/src/_fkset.py',
    'wrapper/src/_lhapdfset.py',
    'wrapper/src/_logger.py',
    'wrapper/src/_nnmpi.py',
    'wrapper/src/_nnpdfdb.py',
    'wrapper/src/_parametrisation.py',
    'wrapper/src/_pathlib.py',
    'wrapper/src/_pdfset.py',
    'wrapper/src/_positivity.py',
    'wrapper/src/_randomgenerator.py',
    'wrapper/src/_thpredictions.py',
    'wrapper/src/_timer.py',
    'wrapper/src/_utils.py',
]

install_data(pyfiles, install_dir: python_install_dir)


wrapper_lib = python.extension_module('_nnpdf',
    [wrapper_cxx],
    dependencies: [wrapper_deps, py3_dep],
	include_directories: inc,
    link_with: [lib],
	install: true,
    install_dir: python_install_dir,
)

#This is to fix https://github.com/ninja-build/ninja/issues/1184
swig_python_target = custom_target('swig_python',
    output: ['@BASENAME@.py'],
	input: [swig_source],
	command:['echo'],
	depends: [wrapper_lib],
	install: true,
    install_dir: python_install_dir,

)

