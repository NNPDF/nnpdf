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

# install scripts
install(FILES ${PROJECT_SOURCE_DIR}/nnpdfcpp/scripts/fitmanager
  DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# install data
install(DIRECTORY ${PROJECT_SOURCE_DIR}/nnpdfcpp/data/ DESTINATION ${PROFILE_PREFIX}/data/)

# install results
install(DIRECTORY DESTINATION ${PROFILE_PREFIX}/results/)

file(WRITE ${PROJECT_BINARY_DIR}/nnpdfcpp/README.md "Store here your fit output.")
install(FILES ${PROJECT_BINARY_DIR}/nnpdfcpp/README.md DESTINATION ${PROFILE_PREFIX}/results/)

# install configs
install(FILES ${PROJECT_SOURCE_DIR}/nnpdfcpp/config/plotting.yml
  ${PROJECT_SOURCE_DIR}/nnpdfcpp/config/fiatlux.yml
  ${PROJECT_SOURCE_DIR}/nnpdfcpp/config/NNPDF30_nnlo_as_0118.yml
  DESTINATION ${PROFILE_PREFIX}/config PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
  )
