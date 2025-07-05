#pragma once

#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/BufferManager.h>
#include <gpu/WebGPU//ComputePipelines.h>
#include <data/atoms/AtomFF.h>

#include <vector>
#include <iostream>
#include <thread>

using namespace ausaxs;
using namespace ausaxs::data;

namespace ausaxs::gpu {
    class WebGPU {
        public:
            WebGPU() {
                initialize();
            }

            void submit_self(const hist::detail::CompactCoordinates& atoms);
            void submit_cross(const hist::detail::CompactCoordinates& atom_1, const hist::detail::CompactCoordinates& atom_2);

            std::vector<double> run();

        private:
            wgpu::Instance instance;
            wgpu::Device device;
            wgpu::Adapter adapter;
            ComputePipelines::Pipelines pipelines;
            BufferManager buffers;
            std::vector<wgpu::Buffer> readback_buffers;

            void initialize();
            wgpu::Instance create_instance();
            wgpu::Adapter get_adapter(wgpu::Instance instance);
            wgpu::Device get_device(wgpu::Instance instance);
    };

    inline void WebGPU::initialize() {
        instance = create_instance();
        device = get_device(instance);
        buffers = BufferManager(device);
        pipelines = ComputePipelines::create(device);
    }

    inline std::vector<double> WebGPU::run() {
        return buffers.merge_histograms(instance);
    }

    inline void WebGPU::submit_self(const hist::detail::CompactCoordinates& atoms) {
        wgpu::CommandEncoder encoder = device.createCommandEncoder();
        wgpu::BindGroup bind_group = buffers.create(encoder, atoms);

        wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
        auto compute_pass = encoder.beginComputePass(compute_pass_desc);
        compute_pass.setPipeline(pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        wgpu::CommandBufferDescriptor command_buffer_desc;
        command_buffer_desc.nextInChain = nullptr;
        auto command_buffer = encoder.finish(command_buffer_desc);
        device.getQueue().submit(command_buffer);
        command_buffer.release();
        encoder.release();
    }

    inline void WebGPU::submit_cross(const hist::detail::CompactCoordinates& atoms1, const hist::detail::CompactCoordinates& atoms2) {
        wgpu::CommandEncoder encoder = device.createCommandEncoder();
        wgpu::BindGroup bind_group = buffers.create(encoder, atoms1, atoms2);

        wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
        auto compute_pass = encoder.beginComputePass(compute_pass_desc);
        compute_pass.setPipeline(pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        wgpu::CommandBufferDescriptor command_buffer_desc;
        command_buffer_desc.nextInChain = nullptr;
        auto command_buffer = encoder.finish(command_buffer_desc);
        device.getQueue().submit(command_buffer);
        command_buffer.release();
        encoder.release();
    }
}