# Include files (should this information not be known at this point?)
include_directories(${PROJECT_SOURCE_DIR}/nnpdfcpp/src/common/inc)
include_directories(${PROJECT_SOURCE_DIR}/nnpdfcpp/src/nnfit/inc)
include_directories(${PROJECT_SOURCE_DIR}/n3fit/evolven3fit)
include_directories(${PROJECT_SOURCE_DIR}/libnnpdf/src/)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/binaries)

configure_file(
  "${PROJECT_SOURCE_DIR}/libnnpdf/src/NNPDF/common.h.in"
  "${PROJECT_SOURCE_DIR}/libnnpdf/src/NNPDF/common.h"
  )

# Add files to the make
add_executable(evolven3fit ${PROJECT_SOURCE_DIR}/n3fit/evolven3fit/evolven3fit.cc
        ${PROJECT_SOURCE_DIR}/nnpdfcpp/src/common/src/md5.cc
        ${PROJECT_SOURCE_DIR}/nnpdfcpp/src/common/src/exportgrid.cc
        ${PROJECT_SOURCE_DIR}/nnpdfcpp/src/nnfit/src/evolgrid.cc )

# Set all flags
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${NNPDF_LDFLAGS} ${GSL_LDFLAGS} ${APFEL_LIBRARIES} ${YAML_LDFLAGS}")

# I am pretty sure this should not be a thing
string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
target_link_libraries(evolven3fit nnpdf ${YAML_LDFLAGS} ${APFEL_LIBRARIES} ${GSL_LDFLAGS})

install(TARGETS evolven3fit DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
