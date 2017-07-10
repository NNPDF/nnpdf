add_executable(catch_test EXCLUDE_FROM_ALL main.cc experiment_test.cc)
target_link_libraries(catch_test nnpdf ${LibArchive_LIBRARIES})
