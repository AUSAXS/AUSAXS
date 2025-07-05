#pragma once

#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/InstanceManager.h>
#include <hist/detail/HistDetailFwd.h>

namespace ausaxs::gpu {
    struct Buffers {
        struct BufferInstance {
            wgpu::Buffer atomic_1;
            wgpu::Buffer atomic_2;
            wgpu::Buffer histogram;
        };

        static BufferInstance create(wgpu::Device device, const hist::detail::CompactCoordinates& atoms);
        static BufferInstance create(wgpu::Device device, const hist::detail::CompactCoordinates& atoms1, const hist::detail::CompactCoordinates& atoms2);
    };
}