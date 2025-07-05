#pragma once

#include <webgpu/webgpu.hpp>

namespace ausaxs::gpu {
    struct ComputePipelines {
        struct Pipelines {
            wgpu::ComputePipeline self;
            wgpu::ComputePipeline cross;
        };

        static ComputePipelines::Pipelines create(wgpu::Device device);
    };
}