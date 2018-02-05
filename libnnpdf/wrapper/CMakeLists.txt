# Set up swig
include(${SWIG_USE_FILE})

include(FindPythonLibs)
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    import numpy as np 
    print(np.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE
  OUTPUT_STRIP_TRAILING_WHITESPACE)

include_directories(${NUMPY_INCLUDE})
include_directories(${PYTHON_INCLUDE_DIRS})

set_source_files_properties(./src/nnpdf.i PROPERTIES CPLUSPLUS ON)
include_directories(../src)

# Add swig module
SWIG_ADD_LIBRARY(nnpdf LANGUAGE python SOURCES ./src/nnpdf.i)
swig_link_libraries(nnpdf ${PROJECT_NAME} ${PYTHON_LIBRARIES})

# Files to install with Python
set(PYTHON_INSTALL_FILES
    ${CMAKE_CURRENT_BINARY_DIR}/nnpdf.py
    ${CMAKE_CURRENT_BINARY_DIR}/_nnpdf.so)

## Configure setup.py and copy to output directory
set(SETUP_PY_IN  ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})

## Install target to call setup.py
add_custom_target(install-python DEPENDS _nnpdf COMMAND python ${SETUP_PY_OUT} install)
