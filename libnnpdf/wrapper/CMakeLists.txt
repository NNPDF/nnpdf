find_package(SWIG REQUIRED)
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs   3 REQUIRED)

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
swig_link_libraries(nnpdf ${PROJECT_NAME} ${PYTHON_LIBRARIES})

install(TARGETS _nnpdf DESTINATION ${PYTHON_SITE_PACKAGES})
install(FILES ${CMAKE_BINARY_DIR}/libnnpdf/wrapper/nnpdf.py DESTINATION ${PYTHON_SITE_PACKAGES})
