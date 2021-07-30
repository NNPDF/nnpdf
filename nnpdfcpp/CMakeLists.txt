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

# install data
install(DIRECTORY ${PROJECT_SOURCE_DIR}/nnpdfcpp/data/ DESTINATION ${PROFILE_PREFIX}/data/ MESSAGE_NEVER)

# install results
install(DIRECTORY DESTINATION ${PROFILE_PREFIX}/results/)
install(DIRECTORY DESTINATION ${PROFILE_PREFIX}/hyperscan_results/)

file(WRITE ${PROJECT_BINARY_DIR}/nnpdfcpp/README.md "Store here your fit output.")
install(FILES ${PROJECT_BINARY_DIR}/nnpdfcpp/README.md DESTINATION ${PROFILE_PREFIX}/results/)

# install configs
install(DIRECTORY ${PROJECT_SOURCE_DIR}/nnpdfcpp/config/ DESTINATION ${PROFILE_PREFIX}/config/)
