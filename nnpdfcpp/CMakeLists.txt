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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# export version to file
configure_file(
  "${PROJECT_SOURCE_DIR}/src/common/inc/version.h.in"
  "${PROJECT_BINARY_DIR}/src/common/inc/version.h"
)
include_directories(${PROJECT_BINARY_DIR}/src/common/inc/)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

find_package(PkgConfig REQUIRED)
pkg_search_module(YAML REQUIRED yaml-cpp)
pkg_search_module(SQLITE3 REQUIRED sqlite3)
pkg_search_module(GSL REQUIRED gsl)
pkg_search_module(NNPDF REQUIRED nnpdf)

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


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NNPDF_CFLAGS} ${LHAPDF_CXX_FLAGS} ${APFEL_CXX_FLAGS} ${YAMLCPP_CXX_FLAGS} ${SQLITE3_CFLAGS} ${GSL_CFLAGS}")

# execute project specific targets
add_subdirectory(src)

# install scripts
install(FILES ${PROJECT_SOURCE_DIR}/scripts/fitmanager
        ${PROJECT_SOURCE_DIR}/scripts/postfit
        DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# install symlin
set(datapath)
execute_process(COMMAND ${PKG_CONFIG_EXECUTABLE} nnpdf --variable=data_path
  OUTPUT_VARIABLE datapath)
string(STRIP "${datapath}" datapath)

if(NOT datapath)
  message(FATAL_ERROR "unexpected: the variable data_path from libnnpdf is empty.")
endif()

if(datapath STREQUAL ${PROJECT_SOURCE_DIR}/data)
  message(WARNING "data path matches")
else()
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/data/ DESTINATION ${datapath})
endif()
