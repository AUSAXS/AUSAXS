#pragma once

#include <gpu/WebGPU/WebGPU.h>

namespace ausaxs::gpu {
    struct ShaderDefinition {
        ShaderDefinition() = default;
        ShaderDefinition(std::string_view source, wgpu::Device device);

        wgpu::ShaderModule module;
        wgpu::BindGroupLayout bind_group_layout;
        struct ComputePipelines {wgpu::ComputePipeline self, cross;} pipelines;
    };
}