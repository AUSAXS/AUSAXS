// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/WebGPU/shaders/ShaderDefinition.h>

// Storage for shaders and their layouts. Compiling and validating a shader is expensive, so each
// definition is created on first use and then shared by all backends using the same device.
namespace ausaxs::gpu::shader {
    struct Simple {
        /// @brief Get the simple distance histogram shader, compiling it on first use.
        template<bool weighted_bins>
        static const ShaderDefinition& get() {
            if constexpr (weighted_bins) {
                return weighted();
            } else {
                return unweighted();
            }
        }

        static const ShaderDefinition& weighted();
        static const ShaderDefinition& unweighted();
    };
}
