cmake_minimum_required (VERSION 3.0.2)

# Disable the use of RPATHS - we probably are not
# that interested in relocatable binaries and it
# causes a lot of trouble on macOS
#SET(CMAKE_SKIP_RPATH TRUE)

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

project(nnpdf)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(VERSION "\"1.2.0b5\"")
set(nnpdfcpp_VERSION 3.1)

# check for dependencies
find_package(PkgConfig REQUIRED)
pkg_search_module(LIBARCHIVE REQUIRED libarchive)
pkg_search_module(SQLITE3 REQUIRED sqlite3)
pkg_search_module(GSL REQUIRED gsl)
pkg_search_module(YAML REQUIRED yaml-cpp)


option(ENABLE_OPENMP "Enable openmp parallelism." OFF)
option(ENABLE_OPENMPI "Enable openMPI parallelism." OFF)
option(ENABLE_TESTS "Enable unit test" OFF)
option(ENABLE_ASAN "Enable ASAN" OFF)
option(VP_DEV "validphys in developer mode" ON)
set(PROFILE_PREFIX "" CACHE STRING "Where you store the 'data' folder. Default empty uses CMAKE_INSTALL_PREFIX/share/NNPDF.")

if (PROFILE_PREFIX)
  set(PROFILE_PREFIX "${PROFILE_PREFIX}")
else (PROFILE_PREFIX)
  set(PROFILE_PREFIX "${CMAKE_INSTALL_PREFIX}/share/NNPDF")
endif()

if (ENABLE_OPENMPI)
  set(LIBNNPDF_HAVE_MPI "#define OPENMPI")
endif()

set(LIBNNPDF_HAVE_SSE "#define SSE_CONV")


if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(ENABLE_OPENMP)

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

set(DEFAULT_CXX_OPTIONS "-Wall -Wextra -march=nocona -mtune=haswell -fvisibility-inlines-hidden -fmessage-length=0 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe")

string(REPLACE "-fsanitize=address" "" CMAKE_EXE_LINKER_FLAGS_STRIPPED "${CMAKE_EXE_LINKER_FLAGS}")
string(REPLACE "-fsanitize=address" "" CMAKE_EXE_LINKER_FLAGS_DEBUG_STRIPPED "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
string(REPLACE "-Wl,-dead_strip_dylibs" "" CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}" CACHE STRING "shared linker flags" FORCE)

# Flags for ASAN
if(ENABLE_ASAN)
  set(CMAKE_CXX_FLAGS "-fsanitize=address ${DEFAULT_CXX_OPTIONS} ${LHAPDF_CXX_FLAGS} ${APFEL_CXX_FLAGS} ${YAML_CFLAGS} ${SQLITE3_CFLAGS} ${GSL_CFLAGS} ${LIBARCHIVE_CFLAGS}" CACHE STRING "compile flags" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g" CACHE STRING "debug compile flags" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS "-fsanitize=address ${CMAKE_EXE_LINKER_FLAGS_STRIPPED}" CACHE STRING "linker flags" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fsanitize=address ${CMAKE_EXE_LINKER_FLAGS_DEBUG_STRIPPED}" CACHE STRING "debug linker flags" FORCE)
else(ENABLE_ASAN)
  set(CMAKE_CXX_FLAGS "${DEFAULT_CXX_OPTIONS} ${LHAPDF_CXX_FLAGS} ${APFEL_CXX_FLAGS} ${YAML_CFLAGS} ${SQLITE3_CFLAGS} ${GSL_CFLAGS} ${LIBARCHIVE_CFLAGS}" CACHE STRING "compile flags" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g" CACHE STRING "debug compile flags" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_STRIPPED}" CACHE STRING "linker flags" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG_STRIPPED}" CACHE STRING "debug linker flags" FORCE)
endif(ENABLE_ASAN)

# libnnpdf configuration
add_subdirectory(libnnpdf)

# nnpdfcpp configuration
add_subdirectory(nnpdfcpp)

# install validphys2
if(VP_DEV)
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install -e ${PROJECT_SOURCE_DIR}/validphys2)")
else(VP_DEV)
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --no-deps --ignore-installed ${PROJECT_SOURCE_DIR}/validphys2)")
endif(VP_DEV)
