cmake_minimum_required (VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# the project name
project(nnpdf)

# activating some global properties for the project
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# export version to file
set(VERSION "\"3.1\"")
configure_file(
  "${PROJECT_SOURCE_DIR}/src/NNPDF/config.h.in"
  "${PROJECT_BINARY_DIR}/src/NNPDF/config.h"
)

option(ENABLE_SAFEMODE "Enable double precision for apfelcomb." OFF)
option(ENABLE_OPENMP "Enable openmp parallelism." OFF)
option(ENABLE_PYWRAP "Enable python wrapper building" ON)

if(NOT ENABLE_SAFEMODE)
  set(LIBNNPDF_HAVE_SSE "#define SSE_CONV")
endif(NOT ENABLE_SAFEMODE)

if(ENABLE_OPENMP)
  set(LIBNNPDF_HAVE_MPI "#define OPENMPI")
endif(ENABLE_OPENMP)

configure_file(
  "${PROJECT_SOURCE_DIR}/wrapper/setup.py.in"
  "${PROJECT_BINARY_DIR}/wrapper/setup.py"
)

set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
configure_file(
  "${PROJECT_SOURCE_DIR}/scripts/nnpdf-config.in"
  "${PROJECT_BINARY_DIR}/scripts/nnpdf-config"
)

# LHAPDF
find_program(LHAPDF_CONFIG lhapdf-config REQUIRED)
if (LHAPDF_CONFIG)
  exec_program(${LHAPDF_CONFIG}
    ARGS --cflags
    OUTPUT_VARIABLE LHAPDF_CXX_FLAGS
  )
  set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${LHAPDF_CONFIG}
    ARGS --libs
    OUTPUT_VARIABLE LHAPDF_LIBRARIES
  )
  set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARIES} CACHE STRING INTERNAL)
endif(LHAPDF_CONFIG)

find_package(GSL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LHAPDF_CXX_FLAGS}")
include_directories(src/NNPDF ${GSL_INCLUDE_DIRS})
add_library(nnpdf SHARED src/common.cc
                         src/commondata.cc
                         src/chisquared.cc
                         src/dataset.cc
                         src/experiments.cc
                         src/fastkernel.cc
                         src/fkconvolute.cc
                         src/fkgenerator.cc
                         src/fkset.cc
                         src/lhapdfset.cc
                         src/logger.cc
                         src/nnmpi.cc
                         src/nnpdfdb.cc
                         src/parametrisation.cc
                         src/pdfset.cc
                         src/positivity.cc
                         src/positivity.cc
                         src/randomgenerator.cc
                         src/thpredictions.cc
                         src/utils.cc
            )

install(FILES ${PROJECT_BINARY_DIR}/scripts/nnpdf-config
	${PROJECT_SOURCE_DIR}/scripts/FKmerge
	${PROJECT_SOURCE_DIR}/scripts/FKselfmerge
	${PROJECT_SOURCE_DIR}/scripts/FKincrement
	DESTINATION bin PERMISSIONS WORLD_READ WORLD_EXECUTE)
install(DIRECTORY src/NNPDF DESTINATION include)
install(TARGETS nnpdf DESTINATION lib)
