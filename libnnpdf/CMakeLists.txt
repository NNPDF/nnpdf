set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")

configure_file(
  "${PROJECT_SOURCE_DIR}/libnnpdf/scripts/nnpdf.pc.in"
  "${PROJECT_SOURCE_DIR}/libnnpdf/scripts/nnpdf.pc"
  )

configure_file(
  "${PROJECT_SOURCE_DIR}/libnnpdf/src/NNPDF/common.h.in"
  "${PROJECT_SOURCE_DIR}/libnnpdf/src/NNPDF/common.h"
  )

# add preprocessor flag
add_definitions(-DDEFAULT_NNPDF_PROFILE_PATH="${prefix}/share/NNPDF/nnprofile.yaml")

# Note that BEFORE is important here: Otherwise we might be reading
# the headers from a previous installation.
include_directories(BEFORE src/NNPDF src)
include_directories(${GSL_INCLUDE_DIRS} ${YAML_INCLUDE_DIRS})
FILE(GLOB_RECURSE Headers "src/NNPDF/*.h")
add_library(nnpdf SHARED src/common.cc
                         src/commondata.cc
                         src/chisquared.cc
                         src/dataset.cc
                         src/experiments.cc
                         src/fastkernel.cc
                         src/fkgenerator.cc
                         src/fkset.cc
                         src/lhapdfset.cc
                         src/logger.cc
                         src/nnmpi.cc
                         src/nnpdfdb.cc
                         src/parametrisation.cc
                         src/pdfset.cc
                         src/positivity.cc
                         src/positivity.cc
                         src/randomgenerator.cc
                         src/thpredictions.cc
                         src/utils.cc
                         src/pathlib.cc
                         ${Headers}
)

target_link_libraries(nnpdf ${LHAPDF_LIBRARIES} ${GSL_LDFLAGS} ${SQLITE3_LDFLAGS} ${LIBARCHIVE_LDFLAGS} ${YAML_LDFLAGS})

install(FILES ${PROJECT_SOURCE_DIR}/libnnpdf/scripts/nnpdf.pc DESTINATION lib/pkgconfig)
install(DIRECTORY src/NNPDF DESTINATION include)
file(WRITE ${PROJECT_SOURCE_DIR}/libnnpdf/REAMDE.md "The share folder of NNPDF has been moved, see: https://github.com/NNPDF/nnpdf/pull/1861\n")
install(FILES ${PROJECT_SOURCE_DIR}/libnnpdf/REAMDE.md DESTINATION share/NNPDF)
install(TARGETS nnpdf DESTINATION lib)
