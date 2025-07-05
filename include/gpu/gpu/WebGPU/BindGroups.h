#pragma once

#include <webgpu/webgpu.hpp>

namespace ausaxs::gpu {
    struct BindGroups {
        static wgpu::BindGroupLayout create(wgpu::Device device);
    };
}