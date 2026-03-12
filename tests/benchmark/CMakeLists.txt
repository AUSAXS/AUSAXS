add_executable(benchmark_pr "benchmark_pr.cpp")
target_link_libraries(benchmark_pr PRIVATE libausaxs Catch2::Catch2WithMain)
target_include_directories(benchmark_pr PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/tests)
set_target_properties(benchmark_pr PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# benchmark histogram calculations using the various simd levels
set(_hist_hot_src ${CMAKE_SOURCE_DIR}/source/core/hist/histogram_manager/HistogramManagerMT.cpp)
set(_hist_includes ${CMAKE_SOURCE_DIR}/include/core ${CMAKE_SOURCE_DIR}/include/math)
set(_base_flags -O3 -ffast-math -fno-finite-math-only)
foreach(_variant scalar sse2 avx2 avx512)
    if(_variant STREQUAL "scalar")
        set(_march -march=native)
        set(_defs   USE_SCALAR)
    elseif(_variant STREQUAL "sse2")
        set(_march -msse2 -mno-avx)
        set(_defs   USE_SSE2)
    elseif(_variant STREQUAL "avx2")
        set(_march -march=x86-64-v3)
        set(_defs   USE_AVX2)
    elseif(_variant STREQUAL "avx512")
        set(_march -march=x86-64-v4)
        set(_defs   USE_AVX512)
    endif()

    # Recompile the hot TU at the target SIMD level
    add_library(hist_hot_${_variant} OBJECT ${_hist_hot_src})
    target_include_directories(hist_hot_${_variant} PRIVATE ${_hist_includes})
    target_compile_options(hist_hot_${_variant} PRIVATE ${_base_flags} ${_march})
    target_compile_definitions(hist_hot_${_variant} PRIVATE ${_defs})

    # Benchmark executable: fresh hot object + full library for the rest
    add_executable(benchmark_histogram_${_variant} "simd_evaluate_benchmark.cpp")
    target_include_directories(benchmark_histogram_${_variant}
        PRIVATE
            ${CMAKE_SOURCE_DIR}/include
            ${CMAKE_SOURCE_DIR}/tests
    )
    target_compile_options(benchmark_histogram_${_variant} PRIVATE ${_base_flags} ${_march})
    target_compile_definitions(benchmark_histogram_${_variant} PRIVATE ${_defs})
    target_link_libraries(benchmark_histogram_${_variant}
        PRIVATE
            $<TARGET_OBJECTS:hist_hot_${_variant}>
            libausaxs_static
            Catch2::Catch2WithMain
    )
    set_target_properties(benchmark_histogram_${_variant} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
    )
endforeach()