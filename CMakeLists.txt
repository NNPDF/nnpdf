cmake_minimum_required (VERSION 3.0.2)

SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "$ENV{CONDA_PREFIX}/lib")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

project(buildmaster)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# check for dependencies
find_package(PkgConfig REQUIRED)
pkg_search_module(GSL REQUIRED gsl)
pkg_search_module(YAML REQUIRED yaml-cpp)
pkg_search_module(NNPDF REQUIRED nnpdf)

set(DEFAULT_CXX_OPTIONS "-Wall -Wextra -march=nocona -mtune=haswell -fvisibility-inlines-hidden -fmessage-length=0 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe")

set(CMAKE_CXX_FLAGS "${DEFAULT_CXX_OPTIONS} ${GSL_CFLAGS} ${YAML_CFLAGS}" CACHE STRING "compile flags" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "linker flags" FORCE)

AUX_SOURCE_DIRECTORY(filters FILTER_SRC)

add_executable(buildmaster src/buildmaster.cc
                       ${FILTER_SRC}
                       src/common.cc
                       src/buildmaster_utils.cc)
  include_directories(src inc)
  target_link_libraries(buildmaster ${NNPDF} ${NNPDF_LDFLAGS} ${YAML_LDFLAGS} ${GSL_LDFLAGS})
  install(TARGETS buildmaster DESTINATION bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)