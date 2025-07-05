#pragma once

#include <webgpu/webgpu.hpp>

namespace ausaxs::gpu {
    struct InstanceManager {
        InstanceManager();
        wgpu::Instance instance;
        wgpu::Device device;

        void process();
    };
}