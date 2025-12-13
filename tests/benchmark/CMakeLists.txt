# Benchmark tests for performance evaluation
# Run with: cmake --build build --target benchmark_tests && ./build/tests/benchmark/bin/benchmark_tests

set(BENCHMARK_TEST_SOURCES
    distance_calculation_benchmark.cpp
)

add_executable(benchmark_tests ${BENCHMARK_TEST_SOURCES} ${CMAKE_SOURCE_DIR}/tests/settings.cpp)

target_link_libraries(benchmark_tests
    PRIVATE
        ausaxs
        Catch2::Catch2WithMain
)

target_include_directories(benchmark_tests
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

set_target_properties(benchmark_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
)

# Don't run benchmarks as part of regular CTest - they're for manual performance testing
# catch_discover_tests(benchmark_tests)
