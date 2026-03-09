add_executable(benchmark_distance_calculation "distance_calculation_benchmark.cpp")

target_link_libraries(benchmark_distance_calculation
    PRIVATE
        libausaxs
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

# Per-SIMD histogram benchmark variants.
#
# Each variant recompiles the hot histogram-manager TU with a specific -march flag so
# the SIMD dispatch macros (__SSE2__, __AVX__, __AVX512F__) reflect that level.  The
# fresh object is linked first; GNU ld then skips the matching archive member from
# libausaxs_static, so no duplicate-symbol issues arise.
#
# Variants:
#   scalar  — -march=native + AUSAXS_FORCE_SCALAR: lets GCC auto-vectorise the
#             scalar fallback, giving a fair comparison against the hand-written paths
#   sse2    — -msse2  (x86-64 baseline with SSE2 only)
#   avx     — -mavx   (AVX without AVX-512)
#   avx512  — -mavx512f

set(_hist_hot_src
    ${CMAKE_SOURCE_DIR}/source/core/hist/histogram_manager/HistogramManagerMT.cpp
)
set(_hist_includes
    ${CMAKE_SOURCE_DIR}/include/core
    ${CMAKE_SOURCE_DIR}/include/math
)
set(_base_flags -O3 -ffast-math -fno-finite-math-only)

foreach(_variant scalar sse2 avx avx512)
    if(_variant STREQUAL "scalar")
        set(_march -march=native)
        set(_defs   AUSAXS_FORCE_SCALAR)
    elseif(_variant STREQUAL "sse2")
        set(_march -msse2 -mno-avx)
        set(_defs   "")
    elseif(_variant STREQUAL "avx")
        set(_march -mavx)
        set(_defs   "")
    elseif(_variant STREQUAL "avx512")
        set(_march -mavx512f)
        set(_defs   "")
    endif()

    # Recompile the hot TU at the target SIMD level
    add_library(hist_hot_${_variant} OBJECT ${_hist_hot_src})
    target_include_directories(hist_hot_${_variant} PRIVATE ${_hist_includes})
    target_compile_options(hist_hot_${_variant} PRIVATE ${_base_flags} ${_march})
    if(_defs)
        target_compile_definitions(hist_hot_${_variant} PRIVATE ${_defs})
    endif()

    # Benchmark executable: fresh hot object + full library for the rest
    add_executable(benchmark_histogram_${_variant} "simd_evaluate_benchmark.cpp")
    target_include_directories(benchmark_histogram_${_variant}
        PRIVATE
            ${CMAKE_SOURCE_DIR}/include
            ${CMAKE_SOURCE_DIR}/tests
    )
    target_compile_options(benchmark_histogram_${_variant} PRIVATE ${_base_flags} ${_march})
    if(_defs)
        target_compile_definitions(benchmark_histogram_${_variant} PRIVATE ${_defs})
    endif()
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
