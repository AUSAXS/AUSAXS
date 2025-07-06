#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/InstanceManager.h>
#include <gpu/WebGPU/BindGroups.h>
#include <gpu/WebGPU/Buffers.h>
#include <gpu/WebGPU/BufferManager.h>
#include <gpu/WebGPU/shaders/ShaderStorage.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/SimpleKernel.h>

#include <vector>
#include <iostream>

using namespace ausaxs;
using namespace ausaxs::data;

namespace ausaxs::gpu {
    template<bool weighted_bins>
    class WebGPU {
        public:
            WebGPU() {
                initialize();
            }

            int submit_self(const hist::detail::CompactCoordinates& atoms, int merge_id = -1);
            int submit_cross(const hist::detail::CompactCoordinates& atom_1, const hist::detail::CompactCoordinates& atom_2, int merge_id = -1);
            hist::distance_calculator::SimpleKernel<weighted_bins>::run_result run();

        private:
            inline static InstanceManager instance;
            Buffers<weighted_bins>::BufferInstance buffers;
            BufferManager<weighted_bins> buffer_manager;
            shader::Simple shaders;

            void initialize();
            wgpu::BindGroup assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer);
    };

    template<bool weighted_bins>
    inline void WebGPU<weighted_bins>::initialize() {
        shaders = shader::Simple(instance.device);
    }

    template<bool weighted_bins>
    inline hist::distance_calculator::SimpleKernel<weighted_bins>::run_result WebGPU<weighted_bins>::run() {
        return buffer_manager.merge(instance);
    }

    template<bool weighted_bins>
    inline wgpu::BindGroup WebGPU<weighted_bins>::assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer) {
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
        bind_group_desc.layout = shaders.get<weighted_bins>().bind_group_layout;
        bind_group_desc.entryCount = entries.size();
        bind_group_desc.entries = entries.data();
        return device.createBindGroup(bind_group_desc);
    }

    template<bool weighted_bins>
    inline int WebGPU<weighted_bins>::submit_self(const hist::detail::CompactCoordinates& atoms, int merge_id) {
        buffers = Buffers<weighted_bins>::create(instance.device, atoms, atoms); //! use single buffer version
        int index = buffer_manager.manage_self(buffers.histogram, merge_id);
        wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atomic_1, buffers.atomic_2, buffers.histogram);

        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();
        auto compute_pass = encoder.beginComputePass(wgpu::Default);
        compute_pass.setPipeline(shaders.get<weighted_bins>().pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        auto command_buffer = encoder.finish();
        instance.device.getQueue().submit(command_buffer);
        command_buffer.release();
        encoder.release();
        return index;
    }

    template<bool weighted_bins>
    inline int WebGPU<weighted_bins>::submit_cross(const hist::detail::CompactCoordinates& atoms1, const hist::detail::CompactCoordinates& atoms2, int merge_id) {
        buffers = Buffers<weighted_bins>::create(instance.device, atoms1, atoms2);
        int index = buffer_manager.manage_cross(buffers.histogram, merge_id);
        wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atomic_1, buffers.atomic_2, buffers.histogram);

        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();
        auto compute_pass = encoder.beginComputePass(wgpu::Default);
        compute_pass.setPipeline(shaders.get<weighted_bins>().pipelines.cross);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        auto command_buffer = encoder.finish();
        instance.device.getQueue().submit(command_buffer);
        command_buffer.release();
        encoder.release();
        return index;
    }
}