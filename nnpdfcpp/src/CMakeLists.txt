# common files
add_library(common STATIC
  common/src/nnpdfsettings.cc
  common/src/md5.cc
  common/src/loadutils.cc
  common/src/datautils.cc
  common/src/common.cc
)
include_directories(common/inc)

# filter
option(COMPILE_filter "Compile filter program." ON)
if (COMPILE_filter)
  add_executable(filter filter/src/filter.cc
                      filter/src/kincuts.cc)
  include_directories(filter/inc)
  target_link_libraries(filter common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES})
endif(COMPILE_filter)

# nnfit
option(COMPILE_nnfit "Compile nnfit program." ON)
if (COMPILE_nnfit)
  add_executable(nnfit nnfit/src/nnfit.cc
                       nnfit/src/apfelevol.cc
                       nnfit/src/fastaddchi2.cc
                       nnfit/src/fitbases.cc
                       nnfit/src/fitpdfset.cc
                       nnfit/src/minimizer.cc
                       nnfit/src/pdfbasis.cc
                       nnfit/src/stopping.cc)
  include_directories(nnfit/inc)
  target_link_libraries(nnfit common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${APFEL_LIBRARIES})
endif(COMPILE_nnfit)

# chi2check
option(COMPILE_chi2check "Compile chi2check program." ON)
if (COMPILE_chi2check)
  add_executable(chi2check chi2check/src/chi2check.cc)
  target_link_libraries(chi2check common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES})
endif(COMPILE_chi2check)

# validphys
option(COMPILE_validphys "Compile validphys program." ON)
if (COMPILE_validphys)
  add_executable(validphys validphys/src/validphys.cc
                           validphys/src/plotdata.cc
                           validphys/src/plotutils.cc)
  include_directories(validphys/inc)
  target_link_libraries(validphys common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${ROOT_LIBRARIES} ${YAMLCPP_LIBRARIES})
endif(COMPILE_validphys)

# fkcheck
option(COMPILE_fkcheck "Compile fkcheck program." OFF)
if (COMPILE_fkcheck)
  add_executable(fkcheck fkcheck/src/fkcheck.cc)
  target_link_libraries(fkcheck common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES})
endif(COMPILE_fkcheck)

# gencovmat
option(COMPILE_gencovmat "Compile gencovmat program." OFF)
if (COMPILE_gencovmat)
  add_executable(gencovmat gencovmat/src/gencovmat.cc)
  target_link_libraries(gencovmat common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES})
endif(COMPILE_gencovmat)

# kinplot
option(COMPILE_kinplot "Compile kinplot program." OFF)
if (COMPILE_kinplot)
  add_executable(kinplot kinplot/src/kinplot.cc)
  target_link_libraries(kinplot common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${ROOT_LIBRARIES})
endif(COMPILE_kinplot)

# revolve
option(COMPILE_revolve "Compile revolve program." OFF)
if (COMPILE_revolve)
  add_executable(revolve revolve/src/revolve.cc)
  target_link_libraries(revolve common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${ROOT_LIBRARIES} ${APFEL_LIBRARIES})
endif(COMPILE_revolve)

# revolve
option(COMPILE_tconvcheck "Compile tconvcheck program." OFF)
if (COMPILE_tconvcheck)
  add_executable(tconvcheck tconvcheck/src/tconvcheck.cc)
  target_link_libraries(tconvcheck common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${ROOT_LIBRARIES} ${APFEL_LIBRARIES})
endif(COMPILE_tconvcheck)
