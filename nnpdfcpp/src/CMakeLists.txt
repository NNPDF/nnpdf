# common files
add_library(common STATIC
  common/src/nnpdfsettings.cc
  common/src/md5.cc
  common/src/loadutils.cc
  common/src/datautils.cc
  common/src/exportgrid.cc
)
include_directories(common/inc)

# filter
option(COMPILE_filter "Compile filter program." OFF)
if (COMPILE_filter)
  add_executable(filter filter/src/filter.cc
                      filter/src/kincuts.cc)
  include_directories(filter/inc)
  target_link_libraries(filter common nnpdf ${LHAPDF_LIBRARIES} ${YAML_LDFLAGS} ${GSL_LDFLAGS})
  install(TARGETS filter DESTINATION bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
endif(COMPILE_filter)

# nnfit
option(COMPILE_nnfit "Compile nnfit program." OFF)
if (COMPILE_nnfit)
  add_executable(nnfit nnfit/src/nnfit.cc
                       nnfit/src/evolgrid.cc
                       nnfit/src/fastaddchi2.cc
                       nnfit/src/fitbases.cc
                       nnfit/src/fitpdfset.cc
                       nnfit/src/minimizer.cc
                       nnfit/src/pdfbasis.cc
                       nnfit/src/stopping.cc)
  include_directories(nnfit/inc)
  target_link_libraries(nnfit common nnpdf ${LHAPDF_LIBRARIES} ${YAML_LDFLAGS} ${APFEL_LIBRARIES} ${GSL_LDFLAGS})
  install(TARGETS nnfit DESTINATION bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
endif(COMPILE_nnfit)

# chi2check
option(COMPILE_chi2check "Compile chi2check program." OFF)
if (COMPILE_chi2check)
  add_executable(chi2check chi2check/src/chi2check.cc)
  target_link_libraries(chi2check common nnpdf ${LHAPDF_LIBRARIES} ${YAML_LDFLAGS} ${GSL_LDFLAGS})
  install(TARGETS chi2check DESTINATION bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
endif(COMPILE_chi2check)


# evolvefit
option(COMPILE_evolvefit "Compile evolvefit program." OFF)
if (COMPILE_evolvefit)
  add_executable(evolvefit
    evolvefit/src/evolvefit.cc
    nnfit/src/evolgrid.cc)
  include_directories(evolvefit/inc nnfit/inc)
  target_link_libraries(evolvefit common nnpdf ${YAML_LDFLAGS} ${APFEL_LIBRARIES} ${GSL_LDFLAGS})
  install(TARGETS evolvefit DESTINATION bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
endif(COMPILE_evolvefit)

# libfiatlux
find_program(FIATLUX_CONFIG fiatlux-config)
if (FIATLUX_CONFIG)
  exec_program(${FIATLUX_CONFIG}
    ARGS --cppflags
    OUTPUT_VARIABLE FIATLUX_CXX_FLAGS
    )
  exec_program(${FIATLUX_CONFIG}
    ARGS --ldflags
    OUTPUT_VARIABLE FIATLUX_LIBRARIES
    )
  # fiatlux
  option(COMPILE_fiatlux "Compile fiatlux program." OFF)
  if (COMPILE_fiatlux)
    add_executable(fiatlux fiatlux/src/fiatlux.cc)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FIATLUX_CXX_FLAGS}")
    target_link_libraries(fiatlux common nnpdf ${LHAPDF_LIBRARIES} ${YAML_LDFLAGS} ${GSL_LDFLAGS} ${FIATLUX_LIBRARIES} ${APFEL_LIBRARIES})
    install(TARGETS fiatlux DESTINATION bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif(COMPILE_fiatlux)
else(FIATLUX_CONFIG)
  message(WARNING "libfiatlux not found: fiatlux compilation is disabled, install libfiatlux and rerun cmake to enable the fiatlux compilation")
endif(FIATLUX_CONFIG)
