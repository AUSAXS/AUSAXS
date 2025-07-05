#pragma once

#include <webgpu/webgpu.hpp>

namespace ausaxs::gpu {
    struct BindGroups {
        /**
         * @brief Get the default bind group layout for the device. 
         */
        static wgpu::BindGroupLayout get(wgpu::Device device);
    };
}