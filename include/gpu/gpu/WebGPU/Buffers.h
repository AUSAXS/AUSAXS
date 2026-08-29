// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/WebGPU/WebGPU.h>
#include <hist/detail/CompactCoordinates.h>

#include <cstdint>

namespace ausaxs::gpu {
    /**
     * @brief Per-dispatch parameters. Must match the Params struct in the wgsl shaders.
     */
    struct ShaderParameters {
        float inv_width;            // inverse bin width; bin = round(inv_width*distance)
        std::uint32_t scale;        // multiplicative factor applied to every contribution
        std::uint32_t bin_count;    // number of bins in the histogram buffer
        std::uint32_t padding = 0;
    };
    static_assert(sizeof(ShaderParameters) == 16, "ShaderParameters must be 16 bytes to match the shader layout.");

    /// @brief The buffers backing a single dispatch.
    struct BufferInstance {
        wgpu::Buffer atomic_1;
        wgpu::Buffer atomic_2;
        wgpu::Buffer histogram;
        wgpu::Buffer parameters;
    };

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

        /// @brief Upload a set of coordinates to a read-only storage buffer.
        template<bool variable_bin_width>
        static wgpu::Buffer create_atom_buffer(const hist::detail::CompactCoordinates<variable_bin_width>& atoms);

        /// @brief Create a zero-initialized histogram buffer with settings::axes::bin_count bins.
        static wgpu::Buffer create_histogram_buffer();

        /// @brief Create a readback buffer matching the size of a histogram buffer.
        static wgpu::Buffer create_readback_buffer();

        /// @brief Upload the per-dispatch parameters to a uniform buffer.
        static wgpu::Buffer create_parameter_buffer(const ShaderParameters& parameters);

        /**
         * @brief Get the shared placeholder buffer used for the unused second atom binding of
         *        self-correlation dispatches. It is owned by the device and must not be released.
         */
        static wgpu::Buffer dummy_buffer();
    };
}
