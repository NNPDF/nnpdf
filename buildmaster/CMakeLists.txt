cmake_minimum_required (VERSION 3.0.2)

# Set a default build type if none was specified

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

project(buildmaster)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
      "${CMAKE_SOURCE_DIR}" CACHE PATH "buildmaster install location" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# Allow people to set their own rpath on command line
SET(CONDA_PREFIX "$ENV{CONDA_PREFIX}")

if((NOT CMAKE_INSTALL_RPATH) AND (CONDA_PREFIX))
  message(STATUS "Setting rpath to $CONDA_PREFIX/lib.")
  SET(CMAKE_INSTALL_RPATH "$ENV{CONDA_PREFIX}/lib")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# check for dependencies

find_package(PkgConfig REQUIRED)
pkg_search_module(GSL REQUIRED gsl)
pkg_search_module(YAML REQUIRED yaml-cpp)
pkg_search_module(NNPDF REQUIRED nnpdf)

# linux complains if we don't explicitly link these

pkg_search_module(SQLITE REQUIRED sqlite3)
pkg_search_module(LHAPDF REQUIRED lhapdf)
pkg_search_module(LIBARCHIVE REQUIRED libarchive)

set(DEFAULT_CXX_OPTIONS "-Wall -Wextra -march=nocona -mtune=haswell \
                         -fvisibility-inlines-hidden -fmessage-length=0 \
                         -ftree-vectorize -fPIC -fstack-protector-strong \
                         -fsanitize=address \
                         -O2 -pipe")

# Get rid of annoying semi colon in pkg_config output - slightly hacky way
# https://stackoverflow.com/a/35458217

string(REPLACE ";" " " NNPDF_CFLAGS "${NNPDF_CFLAGS}")

set(CMAKE_CXX_FLAGS
    "${DEFAULT_CXX_OPTIONS} ${NNPDF_CFLAGS} ${GSL_CFLAGS} ${YAML_CFLAGS}"
    CACHE STRING "compile flags" FORCE)
set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "linker flags" FORCE)

AUX_SOURCE_DIRECTORY(filters FILTER_SRC)

add_executable(buildmaster src/buildmaster.cc
               ${FILTER_SRC}
               src/common.cc
               src/buildmaster_utils.cc)
include_directories(src inc)
target_link_libraries(buildmaster
                      ${NNPDF_LDFLAGS} ${YAML_LDFLAGS}
                      ${GSL_LDFLAGS} ${SQLITE_LDFLAGS}
                      ${LHAPDF_LDFLAGS} ${LIBARCHIVE_LDFLAGS})
install(TARGETS buildmaster DESTINATION .
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
        GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
