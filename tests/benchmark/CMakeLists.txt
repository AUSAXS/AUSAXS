add_executable(benchmark_pr "benchmark_pr.cpp")
target_link_libraries(benchmark_pr PRIVATE ausaxs Catch2::Catch2WithMain)
target_include_directories(benchmark_pr PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/tests)