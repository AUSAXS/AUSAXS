#pragma once

#include <webgpu/webgpu.hpp>

namespace ausaxs::gpu {
    /**
     * @brief This struct initializes the device, and constructs all "one-per-instance" objects.
     */
    struct InstanceManager {
        InstanceManager();
        wgpu::Instance instance;
        wgpu::Device device;

        void process();
    };
}