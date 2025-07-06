#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/InstanceManager.h>
#include <gpu/WebGPU/ComputePipelines.h>
#include <gpu/WebGPU/BindGroups.h>
#include <gpu/WebGPU/Buffers.h>
#include <gpu/WebGPU/BufferManager.h>
#include <hist/detail/CompactCoordinates.h>

#include <vector>
#include <iostream>

using namespace ausaxs;
using namespace ausaxs::data;

namespace ausaxs::gpu {
    class WebGPU {
        public:
            WebGPU() {
                initialize();
            }

            void submit_self(const std::vector<Atom>& atoms, int merge_id = -1);
            void submit_cross(const std::vector<Atom>& atom_1, const std::vector<Atom>& atom_2, int merge_id = -1);

            void submit();

            std::vector<double> run();

        private:
            InstanceManager instance;
            ComputePipelines::Pipelines pipelines;
            Buffers::BufferInstance buffers;
            BufferManager buffer_manager;

            void initialize();
            wgpu::BindGroup assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer, wgpu::BindGroupLayout bind_group_layout);
    };

    inline void WebGPU::initialize() {
        pipelines = ComputePipelines::create(instance);
    }

    inline std::vector<double> WebGPU::run() {
        auto res = buffer_manager.merge(instance);
        return res.self.at(0);
    }

    inline wgpu::BindGroup WebGPU::assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer, wgpu::BindGroupLayout bind_group_layout) {
        assert(buffer1 && "Buffer 1 is null");
        assert(buffer2 && "Buffer 2 is null");
        assert(histogram_buffer && "Histogram buffer is null");
        std::cout << "Assigning buffers to bind group..." << std::endl;

        std::vector<wgpu::BindGroupEntry> entries(3, wgpu::Default);
        entries[0].binding = 0;
        entries[0].buffer = buffer1;
        entries[0].size = buffer1.getSize();

        entries[1].binding = 1;
        entries[1].buffer = buffer2;
        entries[1].size = buffer2.getSize();

        entries[2].binding = 2;
        entries[2].buffer = histogram_buffer;
        entries[2].size = histogram_buffer.getSize();

        wgpu::BindGroupDescriptor bind_group_desc;
        bind_group_desc.layout = bind_group_layout;
        bind_group_desc.entryCount = entries.size();
        bind_group_desc.entries = entries.data();
        return device.createBindGroup(bind_group_desc);
    }

    inline void WebGPU::submit_self(const std::vector<Atom>& atoms, int merge_id) {
        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();

        buffers = Buffers::create(instance.device, atoms);
        buffer_manager.manage_self(buffers.histogram, merge_id);
        wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atomic_1, buffers.atomic_2, buffers.histogram, instance.bind_group_layout);

        // submit work
        auto compute_pass = encoder.beginComputePass(wgpu::Default);
        compute_pass.setPipeline(pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        auto command_buffer = encoder.finish();
        instance.device.getQueue().submit(command_buffer);
        command_buffer.release();
        encoder.release();
    }

    inline void WebGPU::submit_cross(const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2, int merge_id) {
        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();

        buffers = Buffers::create(instance.device, atoms1, atoms2);
        buffer_manager.manage_cross(buffers.histogram, merge_id);
        wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atomic_1, buffers.atomic_2, buffers.histogram, instance.bind_group_layout);

        // submit work
        auto compute_pass = encoder.beginComputePass(wgpu::Default);
        compute_pass.setPipeline(pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        auto command_buffer = encoder.finish();
        instance.device.getQueue().submit(command_buffer);
        command_buffer.release();
        encoder.release();
    }
}