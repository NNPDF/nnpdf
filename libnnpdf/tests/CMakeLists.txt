add_executable(run_test main.cc
                        experiment_test.cc)
target_link_libraries(run_test nnpdf ${LibArchive_LIBRARIES})
