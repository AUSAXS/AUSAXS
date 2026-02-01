add_executable(benchmark_distance_calculation "distance_calculation_benchmark.cpp")

target_link_libraries(benchmark_distance_calculation
    PRIVATE
        ausaxs
        Catch2::Catch2WithMain
)

target_include_directories(benchmark_distance_calculation
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

set_target_properties(benchmark_distance_calculation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
)