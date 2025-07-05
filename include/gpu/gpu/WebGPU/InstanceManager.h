#pragma once

#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/BindGroups.h>

namespace ausaxs::gpu {
    /**
     * @brief This struct initializes the device, and constructs all "one-per-instance" objects.
     */
    struct InstanceManager {
        InstanceManager();
        wgpu::Instance instance;
        wgpu::Device device;
        wgpu::BindGroupLayout bind_group_layout;

        void process();
    };
}