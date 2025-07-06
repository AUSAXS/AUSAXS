#pragma once

#include <gpu/WebGPU/WebGPU.h>
#include <gpu/WebGPU/GPUInstance.h>
#include <hist/detail/HistDetailFwd.h>

namespace ausaxs::gpu {
    template<bool weighted_bins>
    struct Buffers {
        //? Consider containing the webgpu struct definition inside the type definitions, and write them directly into an embedded wgsl script.
        //? This centralizes the type definitions, and makes it significantly easier to ensure they are always in sync.

        // GPU memory layout of unweighted histogram data
        struct HistogramTypeUnweighted {
            float value;

            // for std::plus<> to function with this type
            HistogramTypeUnweighted operator+(const HistogramTypeUnweighted& other) const {
                return {value + other.value};
            }
        };

        // GPU memory layout of weighted histogram data
        struct HistogramTypeWeighted {
            float value;
            std::uint32_t count;
            float bin_center;

            // for std::plus<> to function with this type
            HistogramTypeWeighted operator+(const HistogramTypeWeighted& other) const {
                return {value + other.value, count + other.count, bin_center + other.bin_center};
            }
        };
        using HistogramType = std::conditional_t<weighted_bins, HistogramTypeWeighted, HistogramTypeUnweighted>;

        // Tuple of buffers for use in the GPU pipeline
        struct BufferInstance {
            wgpu::Buffer atomic_1;
            wgpu::Buffer atomic_2;
            wgpu::Buffer histogram;
        };

        static BufferInstance create(wgpu::Device device, const hist::detail::CompactCoordinates& atoms);
        static BufferInstance create(wgpu::Device device, const hist::detail::CompactCoordinates& atoms1, const hist::detail::CompactCoordinates& atoms2);
    };
}