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
set(VERSION "\"1.1.1\"")
configure_file(
  "${PROJECT_SOURCE_DIR}/src/NNPDF/config.h.in"
  "${PROJECT_BINARY_DIR}/src/NNPDF/config.h"
)

option(ENABLE_SAFEMODE "Enable double precision for apfelcomb." OFF)
option(ENABLE_OPENMP "Enable openmp parallelism." OFF)
option(ENABLE_PYWRAP "Enable python wrapper building" OFF)

if(NOT ENABLE_SAFEMODE)
  set(LIBNNPDF_HAVE_SSE "#define SSE_CONV")
  # Check if SSE/AVX instructions are available on the machine where
  # the project is compiled.
  IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
     EXEC_PROGRAM(cat ARGS "/proc/cpuinfo" OUTPUT_VARIABLE CPUINFO)

     STRING(REGEX REPLACE "^.*(sse2).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "sse2" "${SSE_THERE}" SSE2_TRUE)
     IF (SSE2_TRUE)
        set(SSE2_FOUND true CACHE BOOL "SSE2 available on host")
     ELSE (SSE2_TRUE)
        set(SSE2_FOUND false CACHE BOOL "SSE2 available on host")
     ENDIF (SSE2_TRUE)

     # /proc/cpuinfo apparently omits sse3 :(
     STRING(REGEX REPLACE "^.*[^s](sse3).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "sse3" "${SSE_THERE}" SSE3_TRUE)
     IF (NOT SSE3_TRUE)
        STRING(REGEX REPLACE "^.*(T2300).*$" "\\1" SSE_THERE ${CPUINFO})
        STRING(COMPARE EQUAL "T2300" "${SSE_THERE}" SSE3_TRUE)
     ENDIF (NOT SSE3_TRUE)

     STRING(REGEX REPLACE "^.*(ssse3).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "ssse3" "${SSE_THERE}" SSSE3_TRUE)
     IF (SSE3_TRUE OR SSSE3_TRUE)
        set(SSE3_FOUND true CACHE BOOL "SSE3 available on host")
     ELSE (SSE3_TRUE OR SSSE3_TRUE)
        set(SSE3_FOUND false CACHE BOOL "SSE3 available on host")
     ENDIF (SSE3_TRUE OR SSSE3_TRUE)
     IF (SSSE3_TRUE)
        set(SSSE3_FOUND true CACHE BOOL "SSSE3 available on host")
     ELSE (SSSE3_TRUE)
        set(SSSE3_FOUND false CACHE BOOL "SSSE3 available on host")
     ENDIF (SSSE3_TRUE)

     STRING(REGEX REPLACE "^.*(sse4_1).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "sse4_1" "${SSE_THERE}" SSE41_TRUE)
     IF (SSE41_TRUE)
        set(SSE4_1_FOUND true CACHE BOOL "SSE4.1 available on host")
     ELSE (SSE41_TRUE)
        set(SSE4_1_FOUND false CACHE BOOL "SSE4.1 available on host")
     ENDIF (SSE41_TRUE)

     STRING(REGEX REPLACE "^.*(avx).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "avx" "${SSE_THERE}" AVX_TRUE)
     IF (AVX_TRUE)
        set(AVX_FOUND true CACHE BOOL "AVX available on host")
     ELSE (AVX_TRUE)
        set(AVX_FOUND false CACHE BOOL "AVX available on host")
     ENDIF (AVX_TRUE)

     STRING(REGEX REPLACE "^.*(avx2).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "avx2" "${SSE_THERE}" AVX2_TRUE)
     IF (AVX2_TRUE)
        set(AVX2_FOUND true CACHE BOOL "AVX2 available on host")
     ELSE (AVX2_TRUE)
        set(AVX2_FOUND false CACHE BOOL "AVX2 available on host")
     ENDIF (AVX2_TRUE)

  ELSEIF(CMAKE_SYSTEM_NAME MATCHES "Darwin")
     EXEC_PROGRAM("/usr/sbin/sysctl -n machdep.cpu.features" OUTPUT_VARIABLE
        CPUINFO)

     STRING(REGEX REPLACE "^.*[^S](SSE2).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "SSE2" "${SSE_THERE}" SSE2_TRUE)
     IF (SSE2_TRUE)
        set(SSE2_FOUND true CACHE BOOL "SSE2 available on host")
     ELSE (SSE2_TRUE)
        set(SSE2_FOUND false CACHE BOOL "SSE2 available on host")
     ENDIF (SSE2_TRUE)

     STRING(REGEX REPLACE "^.*[^S](SSE3).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "SSE3" "${SSE_THERE}" SSE3_TRUE)
     IF (SSE3_TRUE)
        set(SSE3_FOUND true CACHE BOOL "SSE3 available on host")
     ELSE (SSE3_TRUE)
        set(SSE3_FOUND false CACHE BOOL "SSE3 available on host")
     ENDIF (SSE3_TRUE)

     STRING(REGEX REPLACE "^.*(SSSE3).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "SSSE3" "${SSE_THERE}" SSSE3_TRUE)
     IF (SSSE3_TRUE)
        set(SSSE3_FOUND true CACHE BOOL "SSSE3 available on host")
     ELSE (SSSE3_TRUE)
        set(SSSE3_FOUND false CACHE BOOL "SSSE3 available on host")
     ENDIF (SSSE3_TRUE)

     STRING(REGEX REPLACE "^.*(SSE4.1).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "SSE4.1" "${SSE_THERE}" SSE41_TRUE)
     IF (SSE41_TRUE)
        set(SSE4_1_FOUND true CACHE BOOL "SSE4.1 available on host")
     ELSE (SSE41_TRUE)
        set(SSE4_1_FOUND false CACHE BOOL "SSE4.1 available on host")
     ENDIF (SSE41_TRUE)

     STRING(REGEX REPLACE "^.*(AVX).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "AVX" "${SSE_THERE}" AVX_TRUE)
     IF (AVX_TRUE)
        set(AVX_FOUND true CACHE BOOL "AVX available on host")
     ELSE (AVX_TRUE)
        set(AVX_FOUND false CACHE BOOL "AVX available on host")
     ENDIF (AVX_TRUE)

     STRING(REGEX REPLACE "^.*(AVX2).*$" "\\1" SSE_THERE ${CPUINFO})
     STRING(COMPARE EQUAL "AVX2" "${SSE_THERE}" AVX2_TRUE)
     IF (AVX2_TRUE)
        set(AVX2_FOUND true CACHE BOOL "AVX2 available on host")
     ELSE (AVX2_TRUE)
        set(AVX2_FOUND false CACHE BOOL "AVX2 available on host")
     ENDIF (AVX2_TRUE)

  ELSEIF(CMAKE_SYSTEM_NAME MATCHES "Windows")
     # TODO
     set(SSE2_FOUND   true  CACHE BOOL "SSE2 available on host")
     set(SSE3_FOUND   false CACHE BOOL "SSE3 available on host")
     set(SSSE3_FOUND  false CACHE BOOL "SSSE3 available on host")
     set(SSE4_1_FOUND false CACHE BOOL "SSE4.1 available on host")
     set(AVX_FOUND false CACHE BOOL "AVX available on host")
     set(AVX2_FOUND false CACHE BOOL "AVX2 available on host")
  ELSE(CMAKE_SYSTEM_NAME MATCHES "Linux")
     set(SSE2_FOUND   true  CACHE BOOL "SSE2 available on host")
     set(SSE3_FOUND   false CACHE BOOL "SSE3 available on host")
     set(SSSE3_FOUND  false CACHE BOOL "SSSE3 available on host")
     set(SSE4_1_FOUND false CACHE BOOL "SSE4.1 available on host")
     set(AVX_FOUND false CACHE BOOL "AVX available on host")
     set(AVX2_FOUND false CACHE BOOL "AVX2 available on host")
  ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux")

  if(NOT SSE2_FOUND)
        MESSAGE(STATUS "Could not find hardware support for SSE2 on this machine.")
  endif(NOT SSE2_FOUND)
  if(NOT SSE3_FOUND)
        MESSAGE(STATUS "Could not find hardware support for SSE3 on this machine.")
  endif(NOT SSE3_FOUND)
  if(NOT SSSE3_FOUND)
        MESSAGE(STATUS "Could not find hardware support for SSSE3 on this machine.")
  endif(NOT SSSE3_FOUND)
  if(NOT SSE4_1_FOUND)
        MESSAGE(STATUS "Could not find hardware support for SSE4.1 on this machine.")
  endif(NOT SSE4_1_FOUND)
  if(NOT AVX_FOUND)
        MESSAGE(STATUS "Could not find hardware support for AVX on this machine.")
  endif(NOT AVX_FOUND)
  if(NOT AVX2_FOUND)
        MESSAGE(STATUS "Could not find hardware support for AVX2 on this machine.")
  endif(NOT AVX2_FOUND)

  mark_as_advanced(SSE2_FOUND SSE3_FOUND SSSE3_FOUND SSE4_1_FOUND, AVX_FOUND, AVX2_FOUND)

  if(AVX2_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
  endif()
  if(AVX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
  endif()
  if(SSE4_1_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
  endif()
  if(SSSE3_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mssse3")
  endif()
  if(SSE3_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
  endif()
  if(SSE2_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
  endif()

endif(NOT ENABLE_SAFEMODE)

if(ENABLE_OPENMP)
  set(LIBNNPDF_HAVE_MPI "#define OPENMPI")
  find_package(OpenMP REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(ENABLE_OPENMP)

set(wrappersrc "${PROJECT_SOURCE_DIR}/wrapper/src")
configure_file(
  "${PROJECT_SOURCE_DIR}/wrapper/setup.py.in"
  "${PROJECT_BINARY_DIR}/wrapper/setup.py"
)

set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
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

find_package(GSL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_search_module(SQLITE3 REQUIRED sqlite3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LHAPDF_CXX_FLAGS} ${SQLITE3_CFLAGS}")
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
                         src/utils.cc)
target_link_libraries(nnpdf ${LHAPDF_LIBRARIES} ${GSL_LIBRARIES} ${SQLITE3_LDFLAGS})

install(FILES ${PROJECT_BINARY_DIR}/scripts/nnpdf-config
	${PROJECT_SOURCE_DIR}/scripts/FKmerge
	${PROJECT_SOURCE_DIR}/scripts/FKselfmerge
	${PROJECT_SOURCE_DIR}/scripts/FKincrement
	DESTINATION bin PERMISSIONS WORLD_READ WORLD_EXECUTE)
install(DIRECTORY src/NNPDF DESTINATION include)
install(TARGETS nnpdf DESTINATION lib)

if(ENABLE_PYWRAP)
  add_subdirectory(wrapper)
endif()
