#pragma once

#include <webgpu/webgpu.hpp>

namespace ausaxs::gpu {
    struct GPUInstance {
        GPUInstance();
        wgpu::Instance instance;
        wgpu::Device device;

        void process();
    };
}