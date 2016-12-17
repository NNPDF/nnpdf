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
project(nnpdfcpp)

# activating some global properties for the project
set(nnpdfcpp_VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# export version to file
configure_file(
  "${PROJECT_SOURCE_DIR}/src/common/inc/version.h.in"
  "${PROJECT_BINARY_DIR}/src/common/inc/version.h"
)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

find_package(GSL REQUIRED)
find_package(PkgConfig)
find_package(PkgConfig REQUIRED)
pkg_search_module(YAML REQUIRED yaml-cpp)

find_program(NNPDF_CONFIG nnpdf-config REQUIRED)
if (NNPDF_CONFIG)
  exec_program(${NNPDF_CONFIG}
    ARGS --cppflags
    OUTPUT_VARIABLE NNPDF_CXX_FLAGS
  )
  set(NNPDF_CXX_FLAGS ${NNPDF_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${NNPDF_CONFIG}
    ARGS --ldflags
    OUTPUT_VARIABLE NNPDF_LIBRARIES
  )
  set(NNPDF_LIBRARIES ${NNPDF_LIBRARIES} CACHE STRING INTERNAL)
endif(NNPDF_CONFIG)

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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NNPDF_CXX_FLAGS} ${LHAPDF_CXX_FLAGS}")

# execute project specific targets
add_subdirectory(src)
