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

find_package(PkgConfig REQUIRED)
pkg_search_module(YAML REQUIRED yaml-cpp)
set(YAMLCPP_CXX_FLAGS ${YAML_CFLAGS} CACHE STRING INTERNAL)
set(YAMLCPP_LIBRARIES ${YAML_LDFLAGS} CACHE STRING INTERNAL)
pkg_search_module(SQLITE3 REQUIRED sqlite3)

# NNPDF
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

# APFEL
find_program(APFEL_CONFIG apfel-config REQUIRED)
if (APFEL_CONFIG)
  exec_program(${APFEL_CONFIG}
    ARGS --cppflags
    OUTPUT_VARIABLE APFEL_CXX_FLAGS
  )
  set(APFEL_CXX_FLAGS ${APFEL_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${APFEL_CONFIG}
    ARGS --ldflags
    OUTPUT_VARIABLE APFEL_LIBRARIES
  )
  set(APFEL_LIBRARIES ${APFEL_LIBRARIES} CACHE STRING INTERNAL)
endif(APFEL_CONFIG)

# ROOT (not required, just for validphys)
find_program(ROOT_CONFIG root-config)
if (ROOT_CONFIG)
  exec_program(${ROOT_CONFIG}
    ARGS --cflags
    OUTPUT_VARIABLE ROOT_CXX_FLAGS
  )
  set(ROOT_CXX_FLAGS ${ROOT_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${ROOT_CONFIG}
    ARGS --libs
    OUTPUT_VARIABLE ROOT_LIBRARIES
  )
  set(ROOT_LIBRARIES ${ROOT_LIBRARIES} CACHE STRING INTERNAL)
endif(ROOT_CONFIG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NNPDF_CXX_FLAGS} ${LHAPDF_CXX_FLAGS} ${APFEL_CXX_FLAGS} ${ROOT_CXX_FLAGS} ${YAMLCPP_CXX_FLAGS} ${SQLITE3_CFLAGS}")

# execute project specific targets
add_subdirectory(src)
