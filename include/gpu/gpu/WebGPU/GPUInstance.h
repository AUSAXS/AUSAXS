#pragma once

#include <gpu/WebGPU/WebGPU.h>

namespace ausaxs::gpu {
    struct GPUInstance {
        GPUInstance();
        wgpu::Instance instance;
        wgpu::Device device;

        void process();
        void wait(bool& done);
    };
}