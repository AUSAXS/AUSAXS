// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/WebGPU/WebGPU.h>

#include <string_view>

namespace ausaxs::gpu {
    /**
     * @brief A compiled shader module together with the layout and pipelines needed to dispatch it.
     */
    struct ShaderDefinition {
        ShaderDefinition() = default;

        /**
         * @brief Compile the given wgsl source and create the self- and cross-correlation pipelines.
         * @param source The wgsl source code. See shader::source for the embedded shaders.
         */
        explicit ShaderDefinition(std::string_view source);

        wgpu::ShaderModule module;
        wgpu::BindGroupLayout bind_group_layout;
        struct ComputePipelines {wgpu::ComputePipeline self, cross;} pipelines;
    };
}
