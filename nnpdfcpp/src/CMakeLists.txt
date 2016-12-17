# common files
add_library(common STATIC
  common/src/nnpdfsettings.cc
  common/src/md5.cc
  common/src/loadutils.cc
  common/src/datautils.cc
  common/src/common.cc
)

# filter
add_executable(filter filter/src/filter.cc
                      filter/src/kincuts.cc)
include_directories(common/inc filter/inc)
target_link_libraries(filter common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAML_LIBRARIES})
