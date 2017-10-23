# common files
add_library(common STATIC
  common/src/nnpdfsettings.cc
  common/src/md5.cc
  common/src/loadutils.cc
  common/src/datautils.cc
)
include_directories(common/inc)

# filter
option(COMPILE_filter "Compile filter program." ON)
if (COMPILE_filter)
  add_executable(filter filter/src/filter.cc
                      filter/src/kincuts.cc)
  include_directories(filter/inc)
  target_link_libraries(filter common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${GSL_LDFLAGS})
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
  target_link_libraries(nnfit common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${APFEL_LIBRARIES} ${GSL_LDFLAGS})
endif(COMPILE_nnfit)

# chi2check
option(COMPILE_chi2check "Compile chi2check program." ON)
if (COMPILE_chi2check)
  add_executable(chi2check chi2check/src/chi2check.cc)
  target_link_libraries(chi2check common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${GSL_LDFLAGS})
endif(COMPILE_chi2check)

# ROOT
find_program(ROOT_CONFIG root-config)
if (ROOT_CONFIG)
  exec_program(${ROOT_CONFIG}
    ARGS --cflags
    OUTPUT_VARIABLE ROOT_CXX_FLAGS
    )
  exec_program(${ROOT_CONFIG}
    ARGS --libs
    OUTPUT_VARIABLE ROOT_LIBRARIES
    )
  # validphys
  option(COMPILE_validphys "Compile validphys program." ON)
  if (COMPILE_validphys)  
    add_executable(validphys validphys/src/validphys.cc
      validphys/src/plotdata.cc
      validphys/src/plotutils.cc)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
    include_directories(validphys/inc)
    target_link_libraries(validphys common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${ROOT_LIBRARIES} ${YAMLCPP_LIBRARIES} ${GSL_LDFLAGS})
  endif(COMPILE_validphys)
else(ROOT_CONFIG)
  message(WARNING "ROOT not found: validphys compilation is disabled, install ROOT and rerun cmake to enable the validphys compilation")
endif(ROOT_CONFIG)

# fkcheck
option(COMPILE_fkcheck "Compile fkcheck program." OFF)
if (COMPILE_fkcheck)
  add_executable(fkcheck fkcheck/src/fkcheck.cc)
  target_link_libraries(fkcheck common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${GSL_LDFLAGS})
endif(COMPILE_fkcheck)

# gencovmat
option(COMPILE_gencovmat "Compile gencovmat program." OFF)
if (COMPILE_gencovmat)
  add_executable(gencovmat gencovmat/src/gencovmat.cc)
  target_link_libraries(gencovmat common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${GSL_LDFLAGS})
endif(COMPILE_gencovmat)

# kinplot
option(COMPILE_kinplot "Compile kinplot program." OFF)
if (COMPILE_kinplot)
  add_executable(kinplot kinplot/src/kinplot.cc)
  target_link_libraries(kinplot common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${ROOT_LIBRARIES} ${GSL_LDFLAGS} ${GSL_LDFLAGS})
endif(COMPILE_kinplot)

# revolve
option(COMPILE_revolve "Compile revolve program." OFF)
if (COMPILE_revolve)
  add_executable(revolve revolve/src/revolve.cc)
  target_link_libraries(revolve common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${ROOT_LIBRARIES} ${APFEL_LIBRARIES} ${GSL_LDFLAGS})
endif(COMPILE_revolve)

# revolvenet
option(COMPILE_revolvenet "Compile revolvenet program." OFF)
if (COMPILE_revolvenet)
  add_executable(revolvenet revolvenet/src/revolvenet.cc
    nnfit/src/apfelevol.cc revolvenet/src/nnpdf.cc
    nnfit/src/fitbases.cc)
  include_directories(revolvenet/inc nnfit/inc)
  target_link_libraries(revolvenet common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${ROOT_LIBRARIES} ${APFEL_LIBRARIES} ${GSL_LDFLAGS})
endif(COMPILE_revolvenet)

# revolve
option(COMPILE_tconvcheck "Compile tconvcheck program." OFF)
if (COMPILE_tconvcheck)
  add_executable(tconvcheck tconvcheck/src/tconvcheck.cc)
  target_link_libraries(tconvcheck common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${ROOT_LIBRARIES} ${APFEL_LIBRARIES} ${GSL_LDFLAGS})
endif(COMPILE_tconvcheck)

# mkthpredictions
option(COMPILE_mkthpredictions "Compile mkthpredictions program." OFF)
if (COMPILE_mkthpredictions)
  add_executable(mkthpredictions mkthpredictions/src/mkthpredictions.cc)
  include_directories(mkthpredictions/inc)
  target_link_libraries(mkthpredictions common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${APFEL_LIBRARIES} ${GSL_LDFLAGS})
endif(COMPILE_mkthpredictions)


# filter
option(COMPILE_fiatlux "Compile fiatlux program." OFF)
if (COMPILE_fiatlux)

  # libfiatlux
  find_program(FIATLUX_CONFIG fiatlux-config REQUIRED)
  if (FIATLUX_CONFIG)
    exec_program(${FIATLUX_CONFIG}
      ARGS --cppflags
      OUTPUT_VARIABLE FIATLUX_CXX_FLAGS
    )
    exec_program(${FIATLUX_CONFIG}
      ARGS --ldflags
      OUTPUT_VARIABLE FIATLUX_LIBRARIES
    )
  endif(FIATLUX_CONFIG)

  add_executable(fiatlux fiatlux/src/fiatlux.cc)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FIATLUX_CXX_FLAGS}")
  target_link_libraries(fiatlux common ${NNPDF_LIBRARIES} ${LHAPDF_LIBRARIES} ${YAMLCPP_LIBRARIES} ${GSL_LDFLAGS} ${FIATLUX_LIBRARIES} ${APFEL_LIBRARIES})
endif(COMPILE_fiatlux)
