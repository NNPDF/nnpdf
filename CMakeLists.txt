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
option(VP_DEV "validphys in developer mode" ON)
option(N3_DEV "n3fit in developer mode" ON)
option(BURN_TAG "burn down the git tag in vp and n3fit (only local non-dev installation)" ON)
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

# nnpdfcpp configuration
add_subdirectory(nnpdfcpp)

# evolven3fit
add_subdirectory(n3fit/evolven3fit)

if((BURN_TAG) AND (NOT (VP_DEV) OR NOT (N3_DEV)))
    # Find out where is the root of the git repository (if available) and finds the correct tag to burn down
    # in the package. This just sets the variable, the actual burning is done down below by the validphys
    # or n3fit installation and is found in the version.cmake file as it needs to be done post-installation
    set(GIT_DIR "--git-dir=${PROJECT_SOURCE_DIR}/.git")
    # get the current tag (ex: 3.4) and check at the same time whether git is available and the .git folder found
    execute_process(COMMAND git ${GIT_DIR} describe --abbrev=0 --tags OUTPUT_VARIABLE GIT_TAG ERROR_VARIABLE GIT_ERROR)
    if(GIT_ERROR STREQUAL "")
        string(STRIP ${GIT_TAG} GIT_TAG)
        # get the number of revisions since tag happened
        execute_process(COMMAND git ${GIT_DIR} rev-list ${GIT_TAG}..HEAD --count OUTPUT_VARIABLE GIT_REVN)
        string(STRIP ${GIT_REVN} GIT_REVN)
        # get the shortname for the hash
        execute_process(COMMAND git ${GIT_DIR} rev-parse --short HEAD OUTPUT_VARIABLE GIT_HASH)
        string(STRIP ${GIT_HASH} GIT_HASH)
        # find out whether the repository is in a dirty state
        # (this returns ${GIT_TAG}-dev so we have to remove the tag)
        execute_process(COMMAND git ${GIT_DIR} describe --abbrev=0 --tags --dirty=-dev OUTPUT_VARIABLE GIT_DIRTY)
        string(STRIP ${GIT_DIRTY} GIT_DIRTY)
        string(REPLACE "${GIT_TAG}" "" GIT_DIRTY ${GIT_DIRTY})
        # now concatenate everything
        # with the format ex: 3.4.1880+g00b1741c-dev
        set(GIT_VERSION "build_version='${GIT_TAG}.${GIT_REVN}+g${GIT_HASH}${GIT_DIRTY}'")
        # now set the variable at install time
        install(CODE "set(GIT_VERSION \"${GIT_VERSION}\")")
    else(GIT_ERROR STREQUAL "")
        set(BURN_TAG OFF)
    endif(GIT_ERROR STREQUAL "")
endif()

# install validphys and n3fit
if(VP_DEV)
    install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install -e ${PROJECT_SOURCE_DIR}")
else(VP_DEV)
    install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --no-deps --ignore-installed ${PROJECT_SOURCE_DIR})")
    if(BURN_TAG)
        install(CODE "set(LIBRARY \"validphys\")")
        install(SCRIPT ${PROJECT_SOURCE_DIR}/version.cmake)
        install(CODE "set(LIBRARY \"n3fit\")")
        install(SCRIPT ${PROJECT_SOURCE_DIR}/version.cmake)
    endif(BURN_TAG)
endif(VP_DEV)
