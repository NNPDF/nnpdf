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
	profile_conf.set('PROFILE_PREFIX', join_path(['..', 'nnpdfcpp']))
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

inc = include_directories('src/NNPDF')

lib = shared_library(
    'nnpdf',
	sources,
	include_directories: inc,
	dependencies: dependencies,
	install: true,
    install_dir: 'lib/',
)

pkg_mod = import('pkgconfig')
pkg_mod.generate(
    libraries: lib,
	version: meson.project_version(),
	name: 'libnnpdf',
	filebase: 'nnpdf',
	description: 'NNPDF C++ base library.',

)

#TODO wrappers: This requires some work and upstream fixes.
#swig = find_program('swig')
#swig_sources = [
#    'wrapper/src/dataset.i',
#    'wrapper/src/commondata.i',
#]
#
#wrapper_deps = [
#    dependency('python'),
#    dependency('nnpdf'),
#]
#
#swig_gen = generator(
#    swig,
#	output: ['@BASENAME@.cxx'],
#	arguments: ['-c++', '-python', '-py3', '-o', './@OUTPUT@',
#                '-outdir', '.', '-I../src', '-MMD',
#                '@INPUT@'],
#    depfile: '@BASENAME@.d',
#)
#
#foreach swig_source: swig_sources
#    cpp_source = swig_gen.process(swig_source)
#    shared_library(
#       swig_source.split('.').get(-2),
#       cpp_source,
#       dependencies: dependencies + wrapper_deps,
#       install: true,
#       install_dir: 'lib/',
#    )
#endforeach

#sm = import('swig')
#sm.swig_lib(['lib'], ['wrapper/src/commondata.i'])
