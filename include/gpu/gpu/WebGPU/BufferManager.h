#pragma once

#include <webgpu/webgpu.hpp>

#include <hist/detail/HistDetailFwd.h>

namespace ausaxs::gpu {
    class BufferManager {
        public:
            BufferManager() = default;
            BufferManager(wgpu::Device device) : device(device) {}
            wgpu::BindGroup create(wgpu::CommandEncoder encoder, const hist::detail::CompactCoordinates& atoms);
            wgpu::BindGroup create(wgpu::CommandEncoder encoder, const hist::detail::CompactCoordinates& atoms1, const hist::detail::CompactCoordinates& atoms2);
            std::vector<double> merge_histograms(wgpu::Instance instance) const;

        private:
            wgpu::Device device;
            std::vector<wgpu::Buffer> hist_buffers;
    };
}