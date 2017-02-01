
find_package(PythonInterp 3 REQUIRED)
add_custom_target(wrapper COMMAND mkdir -p ${PROJECT_BINARY_DIR}/wrapper/pylibs/NNPDF
  COMMAND cp ${PROJECT_SOURCE_DIR}/wrapper/src/*.py ${PROJECT_BINARY_DIR}/wrapper/pylibs/NNPDF/
  set(ENV{LD_LIBRARY_PATH}, "/home/zah/anaconda3/envs/swtest/lib/:LD_LIBRARY_PATH")
  COMMAND find ${PROJECT_SOURCE_DIR}/wrapper/src -maxdepth 1 -name "*.i" -print0 | xargs -0 -I {} /home/zah/anaconda3/envs/swtest/bin/swig -c++ -outdir ${PROJECT_BINARY_DIR}/wrapper/pylibs/NNPDF -python -py3 -I${PROJECT_SOURCE_DIR}/src {}
  COMMAND mkdir -p ${PROJECT_BINARY_DIR}/wrapper/_sources
  COMMAND mv ${PROJECT_SOURCE_DIR}/wrapper/src/*.cxx ${PROJECT_BINARY_DIR}/wrapper/_sources
  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_BINARY_DIR}/wrapper/setup.py install
)

add_custom_target(wrapper-clean
  COMMAND rm -rf  ${PROJECT_BINARY_DIR}/_sources
  COMMAND rm -rf  ${PROJECT_BINARY_DIR}/wrapper/build
  COMMAND rm -rf  ${PROJECT_BINARY_DIR}/wrapper/dist
  COMMAND rm -rf  ${PROJECT_BINARY_DIR}/wrapper/pylibs
  COMMAND pip uninstall nnpdf --yes
)
