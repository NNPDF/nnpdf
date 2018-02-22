set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")

configure_file(
  "${PROJECT_SOURCE_DIR}/libnnpdf/nnprofile.yaml.in"
  "${PROJECT_SOURCE_DIR}/libnnpdf/nnprofile.yaml"
  )

configure_file(
  "${PROJECT_SOURCE_DIR}/libnnpdf/src/NNPDF/config.h.in"
  "${PROJECT_SOURCE_DIR}/libnnpdf/src/NNPDF/config.h"
  )

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

include_directories(src/NNPDF ${GSL_INCLUDE_DIRS} ${YAML_INCLUDE_DIRS})
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
            			 src/pathlib.cc)

target_link_libraries(nnpdf ${LHAPDF_LIBRARIES} ${GSL_LDFLAGS} ${SQLITE3_LDFLAGS} ${LIBARCHIVE_LDFLAGS} ${YAML_LDFLAGS})

add_executable(FKconvolute src/fkconvolute.cc)
target_link_libraries(FKconvolute nnpdf ${LibArchive_LIBRARIES})

add_executable(FKmerge2 src/fkmerge2.cc)
target_link_libraries(FKmerge2 nnpdf ${LibArchive_LIBRARIES})

install(FILES ${CMAKE_BINARY_DIR}/libnnpdf/FKconvolute
	${CMAKE_BINARY_DIR}/libnnpdf/FKmerge2
	DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(FILES ${PROJECT_SOURCE_DIR}/libnnpdf/scripts/nnpdf.pc DESTINATION lib/pkgconfig)
install(DIRECTORY src/NNPDF DESTINATION include)
install(FILES ${PROJECT_SOURCE_DIR}/libnnpdf/nnprofile.yaml DESTINATION share/NNPDF)
install(TARGETS nnpdf DESTINATION lib)

add_subdirectory(wrapper)

if(ENABLE_TESTS)
  add_subdirectory(tests)
endif()
