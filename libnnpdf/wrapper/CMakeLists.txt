find_package(SWIG REQUIRED)
find_package(PythonInterp 3.6 REQUIRED)
find_package(PythonLibs   3.6 REQUIRED)

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    import numpy as np 
    print(np.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE
  OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils.sysconfig import get_python_lib; 
    print(get_python_lib())"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES 
  OUTPUT_STRIP_TRAILING_WHITESPACE)

include(${SWIG_USE_FILE})
include_directories(${NUMPY_INCLUDE})
include_directories(${PYTHON_INCLUDE_DIRS})

set_source_files_properties(./src/nnpdf.i PROPERTIES CPLUSPLUS ON)
include_directories(../src)

# Add swig module
if (${CMAKE_VERSION} VERSION_LESS "3.8.0")
    swig_add_module(nnpdf python ./src/nnpdf.i)
else()
    swig_add_library(nnpdf LANGUAGE python SOURCES ./src/nnpdf.i)
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#
https://groups.google.com/a/continuum.io/d/msg/anaconda/057P4uNWyCU/Iem6OtjBCQAJ
set_target_properties(_nnpdf PROPERTIES LINK_FLAGS "-undefined
dynamic_lookup")
swig_link_libraries(nnpdf ${PROJECT_NAME})
else()
swig_link_libraries(nnpdf ${PROJECT_NAME} ${PYTHON_LIBRARIES})
endif()

set(LIBNNPDF_WRAPPER_PREFIX ${PYTHON_SITE_PACKAGES}/NNPDF)
install(TARGETS _nnpdf DESTINATION ${LIBNNPDF_WRAPPER_PREFIX})
install(FILES ${CMAKE_BINARY_DIR}/libnnpdf/wrapper/nnpdf.py DESTINATION ${LIBNNPDF_WRAPPER_PREFIX})
install(FILES ${CMAKE_SOURCE_DIR}/libnnpdf/wrapper/src/__init__.py DESTINATION ${LIBNNPDF_WRAPPER_PREFIX})
