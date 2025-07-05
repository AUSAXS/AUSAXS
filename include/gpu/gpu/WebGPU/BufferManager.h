#pragma once

#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/InstanceManager.h>
#include <hist/detail/HistDetailFwd.h>

namespace ausaxs::gpu {
    class BufferManager {
        public:
            BufferManager(InstanceManager& instance) : instance(instance) {}
            wgpu::BindGroup create(wgpu::CommandEncoder encoder, const hist::detail::CompactCoordinates& atoms);
            wgpu::BindGroup create(wgpu::CommandEncoder encoder, const hist::detail::CompactCoordinates& atoms1, const hist::detail::CompactCoordinates& atoms2);
            std::vector<double> merge_histograms(wgpu::Instance instance) const;

        private:
            InstanceManager& instance;
            std::vector<wgpu::Buffer> hist_buffers;
    };
}