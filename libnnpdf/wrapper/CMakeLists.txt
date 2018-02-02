set(OUTPUT "${PROJECT_BINARY_DIR}/libnnpdf/wrapper/build")

add_custom_command(OUTPUT ${OUTPUT} DEPENDS nnpdf
  COMMAND mkdir -p ${PROJECT_BINARY_DIR}/libnnpdf/wrapper/pylibs/NNPDF
  COMMAND cp ${PROJECT_SOURCE_DIR}/libnnpdf/wrapper/src/*.py ${PROJECT_BINARY_DIR}/libnnpdf/wrapper/pylibs/NNPDF/
  COMMAND find ${PROJECT_SOURCE_DIR}/libnnpdf/wrapper/src -maxdepth 1 -name "*.i" -print0 | xargs -0 -I {} ${SWIG} -c++ -outdir ${PROJECT_BINARY_DIR}/libnnpdf/wrapper/pylibs/NNPDF -python -py3 -I${PROJECT_SOURCE_DIR}/libnnpdf/src {}
  COMMAND mkdir -p ${PROJECT_BINARY_DIR}/libnnpdf/wrapper/_sources
  COMMAND mv ${PROJECT_SOURCE_DIR}/libnnpdf/wrapper/src/*.cxx ${PROJECT_BINARY_DIR}/libnnpdf/wrapper/_sources
  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/libnnpdf/wrapper/setup.py build
  )

add_custom_target(python-wrapper-nnpdf ALL DEPENDS ${OUTPUT})

install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/libnnpdf/wrapper/setup.py install WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/libnnpdf/wrapper)")

add_custom_target(wrapper-clean
  COMMAND rm -rf  ${PROJECT_BINARY_DIR}/libnnpdf/_sources
  COMMAND rm -rf  ${PROJECT_BINARY_DIR}/libnnpdf/wrapper/build
  COMMAND rm -rf  ${PROJECT_BINARY_DIR}/libnnpdfwrapper/dist
  COMMAND rm -rf  ${PROJECT_BINARY_DIR}/libnnpdf/wrapper/pylibs
  COMMAND pip uninstall nnpdf --yes
)
