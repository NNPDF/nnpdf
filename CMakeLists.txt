cmake_minimum_required (VERSION 3.0.2)

# Use rpaths for now, previously there were issues with osx
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

set(VERSION "\"4.0\"")
set(nnpdfcpp_VERSION 4.0)

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
option(ENABLE_DEAD_STRIP "Enable use of flag `-dead_strip-dylibs`" OFF)
option(NNPDF_DEV "n3fit and validphys in developer mode" ON)
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
  #This is to stop the LHAPDF specific warnings from spamming the output
  exec_program(${LHAPDF_CONFIG}
    ARGS --incdir
    OUTPUT_VARIABLE LHAPDF_INCLUDES
  )
  set(LHAPDF_INCLUDES ${LHAPDF_INCLUDES} CACHE STRING INTERNAL)
  include_directories(SYSTEM ${LHAPDF_INCLUDES})
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

#strip linker flags to avoid duplication of asan flags
string(REPLACE "-fsanitize=address" "" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
string(REPLACE "-fsanitize=address" "" CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")

set(CMAKE_ALL_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_MODULE_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")

#dead_strip_dylibs can cause issues, warn if we remove it
if(ENABLE_DEAD_STRIP)
  if(NOT ${CMAKE_ALL_LINKER_FLAGS} MATCHES "-Wl,-dead_strip_dylibs")
    message(WARNING "ENABLE_DEAD_STRIP was set true but there were no occurences of the flag `dead_strip_dylibs` found")
  endif(NOT ${CMAKE_ALL_LINKER_FLAGS} MATCHES "-Wl,-dead_strip_dylibs")
else(ENABLE_DEAD_STRIP)
  if(CMAKE_ALL_LINKER_FLAGS MATCHES "-Wl,-dead_strip_dylibs")
    message(WARNING "removing all occurence of flag: `dead_strip_dylibs` (default behaviour). If you are sure this flag is required then ENABLE_DEAD_STRIP")
    string(REPLACE "-Wl,-dead_strip_dylibs" "" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    string(REPLACE "-Wl,-dead_strip_dylibs" "" CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS}")
    string(REPLACE "-Wl,-dead_strip_dylibs" "" CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
    #executable linker flags get set later
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS}" CACHE STRING "Flags used by the linker during the creation of modules during all build types" FORCE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}" CACHE STRING "Flags used by the linker during the creation of shared libraries during all build types" FORCE)
  endif(CMAKE_ALL_LINKER_FLAGS MATCHES "-Wl,-dead_strip_dylibs")
endif(ENABLE_DEAD_STRIP)

# Flags for ASAN
if(ENABLE_ASAN)
  set(CMAKE_CXX_FLAGS "-fsanitize=address ${DEFAULT_CXX_OPTIONS} ${LHAPDF_CXX_FLAGS} ${APFEL_CXX_FLAGS} ${YAML_CFLAGS} ${SQLITE3_CFLAGS} ${GSL_CFLAGS} ${LIBARCHIVE_CFLAGS}" CACHE STRING "compile flags" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g" CACHE STRING "debug compile flags" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS "-fsanitize=address ${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "linker flags" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fsanitize=address ${CMAKE_EXE_LINKER_FLAGS_DEBUG}" CACHE STRING "debug linker flags" FORCE)
else(ENABLE_ASAN)
  set(CMAKE_CXX_FLAGS "${DEFAULT_CXX_OPTIONS} ${LHAPDF_CXX_FLAGS} ${APFEL_CXX_FLAGS} ${YAML_CFLAGS} ${SQLITE3_CFLAGS} ${GSL_CFLAGS} ${LIBARCHIVE_CFLAGS}" CACHE STRING "compile flags" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g" CACHE STRING "debug compile flags" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "linker flags" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}" CACHE STRING "debug linker flags" FORCE)
endif(ENABLE_ASAN)

# libnnpdf configuration
add_subdirectory(libnnpdf)

# evolven3fit
add_subdirectory(evolven3fit)
install(FILES ${PROJECT_SOURCE_DIR}/validphys2/src/validphys/datafiles/theory.db DESTINATION share/NNPDF)

if(NNPDF_DEV)
    install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --no-deps -e ${PROJECT_SOURCE_DIR})")
else(NNPDF_DEV)
    install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --no-deps --ignore-installed ${PROJECT_SOURCE_DIR})")
endif(NNPDF_DEV)
