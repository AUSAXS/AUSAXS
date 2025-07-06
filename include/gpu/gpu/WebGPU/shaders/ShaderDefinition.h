#pragma once

#include <webgpu/webgpu.hpp>

namespace ausaxs::gpu {
    struct ShaderDefinition {
        ShaderDefinition() = default;
        ShaderDefinition(std::string_view source, wgpu::Device device);

        wgpu::ShaderModule module;
        wgpu::BindGroupLayout bind_group_layout;
        struct ComputePipelines {wgpu::ComputePipeline self, cross;} pipelines;
    };
}