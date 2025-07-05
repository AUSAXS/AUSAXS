#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/InstanceManager.h>
// #include <gpu/WebGPU/BufferManager.h>
#include <gpu/WebGPU/ComputePipelines.h>
#include <gpu/WebGPU/BindGroups.h>

#include <vector>
#include <data/atoms/AtomFF.h>

#include <string_view>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <thread>

using namespace ausaxs;
using namespace ausaxs::data;

namespace ausaxs::gpu {
    struct BufferData {
        wgpu::Buffer atom_1;
        wgpu::Buffer atom_2;
        wgpu::Buffer histogram;
    };

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
            BufferData buffers;
            std::vector<wgpu::Buffer> readback_buffers;

            void initialize();
            wgpu::BindGroup assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer, wgpu::BindGroupLayout bind_group_layout);
            void test_fill_buffers(wgpu::Queue queue, wgpu::Buffer buffer1, wgpu::Buffer buffer2);
            wgpu::ShaderModule load_shader_module(const std::filesystem::path& path, wgpu::Device device);
    };

    BufferData create_buffers_test(wgpu::Device device);
    BufferData create_buffers(wgpu::Device device, const std::vector<Atom>& data);
    BufferData create_buffers(wgpu::Device device, const std::vector<Atom>& data1, const std::vector<Atom>& data2);

    inline void WebGPU::initialize() {
        pipelines = ComputePipelines::create(instance);
    }

    inline void WebGPU::submit() {
        std::cout << "Submitting work..." << std::endl;
        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();

        // fill buffers with test data
        buffers = create_buffers_test(instance.device);
        test_fill_buffers(instance.device.getQueue(), buffers.atom_1, buffers.atom_2);
        wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atom_1, buffers.atom_2, buffers.histogram, instance.bind_group_layout);

        // submit work
        wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
        auto compute_pass = encoder.beginComputePass(compute_pass_desc);
        compute_pass.setPipeline(pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();
        bind_group.release();

        // retrieve data
        static std::string readback_buffer_name = "readback buffer";
        wgpu::BufferDescriptor readback_buffer_desc;
        readback_buffer_desc.label.data = readback_buffer_name.data();
        readback_buffer_desc.label.length = readback_buffer_name.size();
        readback_buffer_desc.size = buffers.histogram.getSize();
        readback_buffer_desc.mappedAtCreation = false;
        readback_buffer_desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
        wgpu::Buffer readback_buffer = instance.device.createBuffer(readback_buffer_desc);
        assert(readback_buffer.getSize() == buffers.histogram.getSize() && "Readback buffer size does not match histogram buffer size.");
        encoder.copyBufferToBuffer(buffers.histogram, 0, readback_buffer, 0, readback_buffer.getSize());

        wgpu::CommandBufferDescriptor command_buffer_desc;
        command_buffer_desc.nextInChain = nullptr;
        command_buffer_desc.label.data = "command buffer";
        auto command_buffer = encoder.finish(command_buffer_desc);
        encoder.release();

        std::cout << "Submitting command buffer..." << std::endl;
        instance.device.getQueue().submit(command_buffer);
        command_buffer.release();
        std::cout << "Command buffer submitted." << std::endl;
        readback_buffers.push_back(readback_buffer);
    }

    inline std::vector<float> WebGPU::run() {
        assert(readback_buffers.size() == 1);
        auto readback_buffer = readback_buffers[0];
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

        std::for_each(readback_buffers.begin(), readback_buffers.end(), [](wgpu::Buffer& buffer) {buffer.release();});
        return result;
    }

    inline wgpu::ShaderModule WebGPU::load_shader_module(const std::filesystem::path& path, wgpu::Device device) {
        std::cout << "Loading shader module..." << std::endl;

        std::ifstream file(path);
        if (!file.is_open()) {
            return nullptr;
        }
        file.seekg(0, std::ios::end);
        size_t size = file.tellg();
        std::string source(size, ' ');
        file.seekg(0);
        file.read(source.data(), size);

        wgpu::ShaderSourceWGSL shader_source;
        shader_source.code = wgpu::StringView(source);
        shader_source.chain.next = nullptr;
        shader_source.chain.sType = wgpu::SType::ShaderSourceWGSL;

        wgpu::ShaderModuleDescriptor shader_module_descriptor{};
        shader_module_descriptor.label = wgpu::StringView("gpu_debug");
        shader_module_descriptor.nextInChain = &shader_source.chain;
        return device.createShaderModule(shader_module_descriptor);
    }

    inline BufferData create_buffers_test(wgpu::Device device) {
        std::cout << "Creating buffers..." << std::endl;
        static std::string atom_buffer_1_name = "atom buffer 1";
        wgpu::BufferDescriptor atom_buffer_1_desc;
        atom_buffer_1_desc.label.data = atom_buffer_1_name.data();
        atom_buffer_1_desc.label.length = atom_buffer_1_name.size();
        atom_buffer_1_desc.size = 8 * sizeof(float);
        atom_buffer_1_desc.mappedAtCreation = false;
        atom_buffer_1_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
        wgpu::Buffer atom_buffer_1 = device.createBuffer(atom_buffer_1_desc);

        static std::string atom_buffer_2_name = "atom buffer 2";
        wgpu::BufferDescriptor atom_buffer_2_desc;
        atom_buffer_2_desc.label.data = atom_buffer_2_name.data();
        atom_buffer_2_desc.label.length = atom_buffer_2_name.size();
        atom_buffer_2_desc.size = 8 * sizeof(float);
        atom_buffer_2_desc.mappedAtCreation = false;
        atom_buffer_2_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
        wgpu::Buffer atom_buffer_2 = device.createBuffer(atom_buffer_2_desc);

        static std::string histogram_buffer_name = "histogram buffer";
        wgpu::BufferDescriptor histogram_buffer_desc;
        histogram_buffer_desc.label.data = histogram_buffer_name.data();
        histogram_buffer_desc.label.length = histogram_buffer_name.size();
        histogram_buffer_desc.size = 8 * sizeof(float);
        histogram_buffer_desc.mappedAtCreation = false;
        histogram_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
        wgpu::Buffer histogram_buffer = device.createBuffer(histogram_buffer_desc);
        std::cout << "Buffers created." << std::endl;

        return {atom_buffer_1, atom_buffer_2, histogram_buffer};
    }

    inline std::vector<float> convert_data(const std::vector<Atom>& atoms) {
        std::vector<float> data(atoms.size()*4);
        for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
            const auto& atom = atoms[i];
            data[i*4 + 0] = atom.x();
            data[i*4 + 1] = atom.y();
            data[i*4 + 2] = atom.z();
            data[i*4 + 3] = atom.w;
        }
        return data;
    }

    inline wgpu::Buffer create_buffer(wgpu::Device device, const std::vector<Atom>& atoms) {
        auto queue = device.getQueue();
        auto data1 = convert_data(atoms);
        wgpu::BufferDescriptor atom_buffer_desc;
        atom_buffer_desc.size = data1.size()*sizeof(float);
        atom_buffer_desc.mappedAtCreation = false;
        atom_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
        wgpu::Buffer atom_buffer = device.createBuffer(atom_buffer_desc);
        queue.writeBuffer(atom_buffer, 0, data1.data(), atom_buffer.getSize());
        return atom_buffer;
    }

    inline BufferData create_buffers(wgpu::Device device, const std::vector<Atom>& atoms) {
        assert(false && "Not implemented yet");
        auto atom_buffer_1 = create_buffer(device, atoms);

        wgpu::BufferDescriptor histogram_buffer_desc;
        histogram_buffer_desc.size = constants::axes::d_vals.size()*sizeof(float);
        histogram_buffer_desc.mappedAtCreation = false;
        histogram_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
        wgpu::Buffer histogram_buffer = device.createBuffer(histogram_buffer_desc);
        return {atom_buffer_1, {}, histogram_buffer};
    }

    inline BufferData create_buffers(wgpu::Device device, const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2) {
        auto atom_buffer_1 = create_buffer(device, atoms1);
        auto atom_buffer_2 = create_buffer(device, atoms2);

        wgpu::BufferDescriptor histogram_buffer_desc;
        histogram_buffer_desc.size = constants::axes::d_vals.size()*sizeof(float);
        histogram_buffer_desc.mappedAtCreation = false;
        histogram_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
        wgpu::Buffer histogram_buffer = device.createBuffer(histogram_buffer_desc);
        return {atom_buffer_1, atom_buffer_2, histogram_buffer};
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

    inline void WebGPU::test_fill_buffers(wgpu::Queue queue, wgpu::Buffer buffer1, wgpu::Buffer buffer2) {
        std::cout << "Filling buffers..." << std::endl;
        std::array<float, 8> data1 = {0, 0, 1, 2, 0, 0, 2, 4};
        std::array<float, 8> data2 = {1, 0, 1, 6, 1, 0, 2, 8};
        queue.writeBuffer(buffer1, 0, data1.data(), buffer1.getSize());
        queue.writeBuffer(buffer2, 0, data2.data(), buffer2.getSize());
    }

    inline void WebGPU::submit_self(const std::vector<Atom>& atoms) {
        std::cout << "Calculating self..." << std::endl;
        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();

        buffers = create_buffers(instance.device, atoms, atoms); //! use single buffer version for self calculation
        wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atom_1, buffers.atom_2, buffers.histogram, instance.bind_group_layout);

        // submit work
        wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
        auto compute_pass = encoder.beginComputePass(compute_pass_desc);
        compute_pass.setPipeline(pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        // retrieve data
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
        readback_buffers.push_back(readback_buffer);    
    }

    inline void WebGPU::submit_cross(const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2) {
        std::cout << "Calculating cross..." << std::endl;
        wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();

        buffers = create_buffers(instance.device, atoms1, atoms2);
        wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atom_1, buffers.atom_2, buffers.histogram, instance.bind_group_layout);

        // submit work
        wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
        auto compute_pass = encoder.beginComputePass(compute_pass_desc);
        compute_pass.setPipeline(pipelines.self);
        compute_pass.setBindGroup(0, bind_group, 0, nullptr);
        compute_pass.dispatchWorkgroups(1, 1, 1);
        compute_pass.end();

        // retrieve data
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
        readback_buffers.push_back(readback_buffer);
    }
}