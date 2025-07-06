#pragma once

#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/InstanceManager.h>
#include <gpu/WebGPU/shaders/ShaderDefinition.h>

namespace ausaxs::gpu {
    template<bool weighted_bins>
    struct ComputePipelines {
        struct Pipelines {
            wgpu::ComputePipeline self;
            wgpu::ComputePipeline cross;
        };

        static ComputePipelines::Pipelines create(wgpu::Device device, ShaderDefinition shader);
    };
}