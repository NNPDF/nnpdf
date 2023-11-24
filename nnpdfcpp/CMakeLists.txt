# export version to file
configure_file(
  "${PROJECT_SOURCE_DIR}/nnpdfcpp/src/common/inc/version.h.in"
  "${PROJECT_BINARY_DIR}/nnpdfcpp/src/common/inc/version.h"
)
include_directories(${PROJECT_BINARY_DIR}/nnpdfcpp/src/common/inc/)
include_directories(${PROJECT_SOURCE_DIR}/libnnpdf/src/)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/binaries)

# execute project specific targets
add_subdirectory(src)

# install configs
install(DIRECTORY ${PROJECT_SOURCE_DIR}/nnpdfcpp/config/ DESTINATION ${PROFILE_PREFIX}/config/)
