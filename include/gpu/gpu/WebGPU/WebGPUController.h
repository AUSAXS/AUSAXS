#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/InstanceManager.h>
#include <gpu/WebGPU/ComputePipelines.h>
#include <gpu/WebGPU/BindGroups.h>
#include <gpu/WebGPU/Buffers.h>
#include <hist/detail/CompactCoordinates.h>

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

            void submit_self(const std::vector<Atom>& atoms);
            void submit_cross(const std::vector<Atom>& atom_1, const std::vector<Atom>& atom_2);

            void submit();

            std::vector<float> run();

        private:
            InstanceManager instance;
            ComputePipelines::Pipelines pipelines;
            Buffers::BufferInstance buffers;

            void initialize();
            wgpu::BindGroup assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer, wgpu::BindGroupLayout bind_group_layout);
    };

    inline void WebGPU::initialize() {
        pipelines = ComputePipelines::create(instance);
    }

    inline std::vector<float> WebGPU::run() {
        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();
        wgpu::BufferDescriptor readback_buffer_desc;
        readback_buffer_desc.size = buffers.histogram.getSize();
        readback_buffer_desc.mappedAtCreation = false;
        readback_buffer_desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
        wgpu::Buffer readback_buffer = instance.device.createBuffer(readback_buffer_desc);

        assert(readback_buffer.getSize() == buffers.histogram.getSize() && "Readback buffer size does not match histogram buffer size.");
        encoder.copyBufferToBuffer(buffers.histogram, 0, readback_buffer, 0, readback_buffer.getSize());
        wgpu::CommandBufferDescriptor command_buffer_desc;
        command_buffer_desc.nextInChain = nullptr;
        auto command_buffer = encoder.finish(command_buffer_desc);
        encoder.release();

        instance.device.getQueue().submit(command_buffer);
        command_buffer.release();

        static std::vector<float> result(readback_buffer.getSize());
        bool done = false;
        wgpu::BufferMapCallbackInfo map_callback;
        map_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
        map_callback.userdata1 = &readback_buffer;
        map_callback.userdata2 = &done;
        map_callback.callback = [] (WGPUMapAsyncStatus, WGPUStringView, void* p_readback_buffer, void* p_done) -> void {
            std::cout << "Readback buffer callback!" << std::endl;
            wgpu::Buffer readback_buffer = *reinterpret_cast<wgpu::Buffer*>(p_readback_buffer);
            bool* done = reinterpret_cast<bool*>(p_done);
            const float* output = (const float*) readback_buffer.getConstMappedRange(0, readback_buffer.getSize());
            result.assign(output, output + readback_buffer.getSize());
            *done = true;
            readback_buffer.unmap();
        };
        readback_buffer.mapAsync(wgpu::MapMode::Read, 0, readback_buffer.getSize(), map_callback);
        instance.process();

        while(!done) {
            std::cout << "Waiting for readback..." << std::endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            instance.process();
        }

        return result;
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
        bind_group_desc.entries = (WGPUBindGroupEntry*) entries.data();
        return device.createBindGroup(bind_group_desc);
    }

    inline void WebGPU::submit_self(const std::vector<Atom>& atoms) {
        std::cout << "Calculating self..." << std::endl;
        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();

        buffers = Buffers::create(instance.device, atoms, atoms); //! use single buffer version for self calculation
        wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atomic_1, buffers.atomic_2, buffers.histogram, instance.bind_group_layout);

        // submit work
        wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
        auto compute_pass = encoder.beginComputePass(compute_pass_desc);
        compute_pass.setPipeline(pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        wgpu::CommandBufferDescriptor command_buffer_desc;
        command_buffer_desc.nextInChain = nullptr;
        auto command_buffer = encoder.finish(command_buffer_desc);
        encoder.release();

        instance.device.getQueue().submit(command_buffer);
        command_buffer.release();
    }

    inline void WebGPU::submit_cross(const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2) {
        std::cout << "Calculating cross..." << std::endl;
        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();

        buffers = Buffers::create(instance.device, atoms1, atoms2);
        wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atomic_1, buffers.atomic_2, buffers.histogram, instance.bind_group_layout);

        // submit work
        wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
        auto compute_pass = encoder.beginComputePass(compute_pass_desc);
        compute_pass.setPipeline(pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        wgpu::CommandBufferDescriptor command_buffer_desc;
        command_buffer_desc.nextInChain = nullptr;
        auto command_buffer = encoder.finish(command_buffer_desc);
        encoder.release();

        instance.device.getQueue().submit(command_buffer);
        command_buffer.release();
    }
}