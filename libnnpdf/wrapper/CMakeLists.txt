find_package(SWIG 3.0.10 REQUIRED)
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs   3 REQUIRED)

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    import numpy as np
    print(np.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE
	RESULT_VARIABLE ret
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT ret EQUAL "0")
	message(FATAL_ERROR "numpy not found")
endif()

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils.sysconfig import get_python_lib;
    print(get_python_lib(prefix=\"${CMAKE_INSTALL_PREFIX}\"))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
	RESULT_VARIABLE ret
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT ret EQUAL "0")
	message(FATAL_ERROR "Could not get python compiler info")
endif()

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    import sysconfig;
    lcom = sysconfig.get_config_var(\"LDSHARED\");
    cc, _, flags = lcom.strip().partition(\" \");
    print(flags)"
    OUTPUT_VARIABLE PYTHON_EXTENSION_LINKFLAGS
	RESULT_VARIABLE ret
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT ret EQUAL "0")
	message(FATAL_ERROR "Could not get python sysconfig info")
endif()

include(${SWIG_USE_FILE})
include_directories(${NUMPY_INCLUDE})
include_directories(${PYTHON_INCLUDE_DIRS})

set_source_files_properties(./src/nnpdf.i PROPERTIES CPLUSPLUS ON)
include_directories(../src)

# Add swig module
SET(CMAKE_SWIG_FLAGS "-threads")
if (${CMAKE_VERSION} VERSION_LESS "3.8.0")
    swig_add_module(nnpdf python ./src/nnpdf.i)
else()
    swig_add_library(nnpdf LANGUAGE python SOURCES ./src/nnpdf.i)
endif()
swig_link_libraries(nnpdf ${PROJECT_NAME} ${PYTHON_EXTENSION_LINKFLAGS})

set(LIBNNPDF_WRAPPER_PREFIX ${PYTHON_SITE_PACKAGES}/NNPDF)
install(TARGETS _nnpdf DESTINATION ${LIBNNPDF_WRAPPER_PREFIX})
install(FILES ${CMAKE_BINARY_DIR}/libnnpdf/wrapper/nnpdf.py DESTINATION ${LIBNNPDF_WRAPPER_PREFIX})
install(FILES ${CMAKE_SOURCE_DIR}/libnnpdf/wrapper/src/__init__.py DESTINATION ${LIBNNPDF_WRAPPER_PREFIX})
