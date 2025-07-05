#include "hist/detail/CompactCoordinates.h"
#include <webgpu/webgpu.hpp>

// #include <vector>
// #include <data/atoms/AtomFF.h>

// #include <string_view>
// #include <iostream>
// #include <filesystem>
// #include <fstream>
// #include <thread>

// using namespace ausaxs;
// using namespace ausaxs::data;

// struct BufferData {
//     wgpu::Buffer atom_1;
//     wgpu::Buffer atom_2;
//     wgpu::Buffer histogram;
// };

// struct ComputePipelines {
//     wgpu::ComputePipeline self;
//     wgpu::ComputePipeline cross;
// };

// class WebGPU {
//     public:
//         WebGPU() {
//             initialize();
//         }

//         void submit_self(const std::vector<Atom>& atoms);
//         void submit_cross(const std::vector<Atom>& atom_1, const std::vector<Atom>& atom_2);

//         void submit();

//         std::vector<float> run();

//     private:
//         wgpu::Instance instance;
//         wgpu::Device device;
//         wgpu::Adapter adapter;
//         wgpu::BindGroupLayout bind_group_layout;
//         ComputePipelines pipelines;
//         BufferData buffers;
//         std::vector<wgpu::Buffer> readback_buffers;

//         void initialize();
//         wgpu::Instance create_instance();
//         wgpu::Adapter get_adapter(wgpu::Instance instance);
//         wgpu::Device get_device(wgpu::Instance instance);
//         wgpu::BindGroupLayout initialize_bind_group_layout(wgpu::Device device);
//         ComputePipelines initialize_compute_pipelines(wgpu::Device device, wgpu::BindGroupLayout bind_group_layout);
//         wgpu::BindGroup assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer, wgpu::BindGroupLayout bind_group_layout);
//         void test_fill_buffers(wgpu::Queue queue, wgpu::Buffer buffer1, wgpu::Buffer buffer2);
//         wgpu::ShaderModule load_shader_module(const std::filesystem::path& path, wgpu::Device device);
// };

// BufferData create_buffers_test(wgpu::Device device);
// BufferData create_buffers(wgpu::Device device, const std::vector<Atom>& data);
// BufferData create_buffers(wgpu::Device device, const std::vector<Atom>& data1, const std::vector<Atom>& data2);

// void WebGPU::initialize() {
//     instance = create_instance();
//     device = get_device(instance);
//     bind_group_layout = initialize_bind_group_layout(device);
//     pipelines = initialize_compute_pipelines(device, bind_group_layout);
// }

// void WebGPU::submit() {
//     std::cout << "Submitting work..." << std::endl;
//     wgpu::CommandEncoder encoder = device.createCommandEncoder();

//     // fill buffers with test data
//     buffers = create_buffers_test(device);
//     test_fill_buffers(device.getQueue(), buffers.atom_1, buffers.atom_2);
//     wgpu::BindGroup bind_group = assign_buffers(device, buffers.atom_1, buffers.atom_2, buffers.histogram, bind_group_layout);

//     // submit work
//     wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
//     auto compute_pass = encoder.beginComputePass(compute_pass_desc);
//     compute_pass.setPipeline(pipelines.self);
//     compute_pass.setBindGroup(0, bind_group, 0, nullptr);
//     compute_pass.dispatchWorkgroups(1, 1, 1);
//     compute_pass.end();
//     bind_group.release();

//     // retrieve data
//     static std::string readback_buffer_name = "readback buffer";
//     wgpu::BufferDescriptor readback_buffer_desc;
//     readback_buffer_desc.label.data = readback_buffer_name.data();
//     readback_buffer_desc.label.length = readback_buffer_name.size();
//     readback_buffer_desc.size = buffers.histogram.getSize();
//     readback_buffer_desc.mappedAtCreation = false;
//     readback_buffer_desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
//     wgpu::Buffer readback_buffer = device.createBuffer(readback_buffer_desc);
//     assert(readback_buffer.getSize() == buffers.histogram.getSize() && "Readback buffer size does not match histogram buffer size.");
//     encoder.copyBufferToBuffer(buffers.histogram, 0, readback_buffer, 0, readback_buffer.getSize());

//     wgpu::CommandBufferDescriptor command_buffer_desc;
//     command_buffer_desc.nextInChain = nullptr;
//     command_buffer_desc.label.data = "command buffer";
//     auto command_buffer = encoder.finish(command_buffer_desc);
//     encoder.release();

//     std::cout << "Submitting command buffer..." << std::endl;
//     device.getQueue().submit(command_buffer);
//     command_buffer.release();
//     std::cout << "Command buffer submitted." << std::endl;
//     readback_buffers.push_back(readback_buffer);
// }

// std::vector<float> WebGPU::run() {
//     assert(readback_buffers.size() == 1);
//     auto readback_buffer = readback_buffers[0];
//     static std::vector<float> result(readback_buffer.getSize());

//     bool done = false;
//     wgpu::BufferMapCallbackInfo map_callback;
//     map_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
//     map_callback.userdata1 = &readback_buffer;
//     map_callback.userdata2 = &done;
//     map_callback.callback = [] (WGPUMapAsyncStatus, WGPUStringView, void* p_readback_buffer, void* p_done) -> void {
//         std::cout << "Readback buffer callback!" << std::endl;
//         wgpu::Buffer readback_buffer = *reinterpret_cast<wgpu::Buffer*>(p_readback_buffer);
//         bool* done = reinterpret_cast<bool*>(p_done);
//         const float* output = (const float*) readback_buffer.getConstMappedRange(0, readback_buffer.getSize());
//         result.assign(output, output + readback_buffer.getSize());
//         *done = true;
//         readback_buffer.unmap();
//     };
//     readback_buffer.mapAsync(wgpu::MapMode::Read, 0, readback_buffer.getSize(), map_callback);
//     instance.processEvents();

//     while(!done) {
//         std::cout << "Waiting for readback..." << std::endl;
//         std::this_thread::sleep_for(std::chrono::milliseconds(100));
//         instance.processEvents();
//     }

//     std::for_each(readback_buffers.begin(), readback_buffers.end(), [](wgpu::Buffer& buffer) {buffer.release();});
//     return result;
// }

// wgpu::ShaderModule WebGPU::load_shader_module(const std::filesystem::path& path, wgpu::Device device) {
//     std::cout << "Loading shader module..." << std::endl;

//     std::ifstream file(path);
//     if (!file.is_open()) {
//         return nullptr;
//     }
//     file.seekg(0, std::ios::end);
//     size_t size = file.tellg();
//     std::string source(size, ' ');
//     file.seekg(0);
//     file.read(source.data(), size);

//     wgpu::ShaderSourceWGSL shader_source;
//     shader_source.code = wgpu::StringView(source);
//     shader_source.chain.next = nullptr;
//     shader_source.chain.sType = wgpu::SType::ShaderSourceWGSL;

//     wgpu::ShaderModuleDescriptor shader_module_descriptor{};
//     shader_module_descriptor.label = wgpu::StringView("gpu_debug");
//     shader_module_descriptor.nextInChain = &shader_source.chain;
//     return device.createShaderModule(shader_module_descriptor);
// }

// wgpu::Instance WebGPU::create_instance() {
//     std::cout << "Creating WebGPU instance..." << std::endl;
//     wgpu::InstanceDescriptor descriptor;
//     wgpu::Instance instance = wgpu::createInstance(descriptor);
//     assert(instance && "Could not create WebGPU instance.");
//     return instance;
// }

// wgpu::Adapter WebGPU::get_adapter(wgpu::Instance instance) {
//     if (!instance) {
//         std::cout << "Could not create WebGPU instance." << std::endl;
//         exit(0);
//     } else {
//         std::cout << "WebGPU instance created." << std::endl;
//     }

//     struct UserData {
//         wgpu::Adapter adapter;
//         bool request_finished = false;
//     }; 
//     UserData user_data;

//     {   // perform the adapter request        
//         wgpu::RequestAdapterOptions options;
//         options.powerPreference = wgpu::PowerPreference::HighPerformance;

//         // define callback function triggered when the adapter is ready (or failed)
//         wgpu::RequestAdapterCallbackInfo adapter_callback;
//         adapter_callback.userdata1 = &user_data;
//         adapter_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
//         adapter_callback.callback = [] (WGPURequestAdapterStatus status, WGPUAdapter adapter, WGPUStringView message, void* user_data_p, void*) -> void {
//             UserData& user_data = *reinterpret_cast<UserData*>(user_data_p);
//             if (status == wgpu::RequestAdapterStatus::Success) {
//                 std::cout << "Adapter found." << std::endl;
//                 user_data.adapter = adapter;
//             } else {
//                 std::cout << "Could not get WebGPU adapter: " << std::string(message.data, message.length) << std::endl;
//             }
//             user_data.request_finished = true;
//         };
//         instance.requestAdapter(options, adapter_callback);
//         instance.processEvents();
//     }

//     assert(user_data.request_finished);

//     // print adapter info
//     wgpu::AdapterInfo info;
//     user_data.adapter.getInfo(&info);
//     std::cout << "WebGPU: Acquired device: " << std::string(info.device.data, info.device.length) << std::endl;

//     return user_data.adapter;
// }

// wgpu::Device WebGPU::get_device(wgpu::Instance instance) {
//     // create a device with default descriptor
//     auto adapter = get_adapter(instance);
//     wgpu::DeviceDescriptor device_descriptor;

//     {   
//         // setup callback function for disconnected devices
//         wgpu::DeviceLostCallbackInfo device_callback;
//         device_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
//         device_callback.callback = [] (WGPUDeviceImpl *const *, WGPUDeviceLostReason reason, WGPUStringView message, void *, void *) -> void {
//             std::cout << "Device lost: reason " << reason;
//             if (message.data) std::cout << " (" << std::string(message.data, message.length) << ")";
//             std::cout << std::endl;
//         };

//         // setup callback function for uncaptured errors
//         wgpu::UncapturedErrorCallbackInfo error_callback;
//         error_callback.callback = [] (WGPUDeviceImpl *const *, WGPUErrorType error, WGPUStringView message, void *, void *) -> void {
//             std::cout << "Uncaptured error: " << error;
//             if (message.data) std::cout << " (" << std::string(message.data, message.length) << ")";
//             std::cout << std::endl;
//         };

//         device_descriptor.uncapturedErrorCallbackInfo = error_callback;
//         device_descriptor.deviceLostCallbackInfo = device_callback;
//     }

//     struct UserData {
//         WGPUDevice device = nullptr;
//         bool requestEnded = false;
//     };
//     UserData userData;

//     {   // perform the device request
//         WGPURequestDeviceCallbackInfo device_callback;
//         device_callback.userdata1 = &userData;
//         device_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
//         device_callback.callback = [] (WGPURequestDeviceStatus status, WGPUDevice device, WGPUStringView message, void* user_data_p, void*) -> void {
//             UserData& userData = *reinterpret_cast<UserData*>(user_data_p);
//             if (status == wgpu::RequestDeviceStatus::Success) {
//                 std::cout << "Device found." << std::endl;
//                 userData.device = device;
//             } else {
//                 std::cout << "Could not get WebGPU device: " << std::string(message.data, message.length) << std::endl;
//             }
//             userData.requestEnded = true;
//         };
//         adapter.requestDevice(device_descriptor, device_callback);
//         instance.processEvents();
//     }

//     assert(userData.requestEnded);
//     adapter.release();
//     return userData.device;
// }

// wgpu::BindGroupLayout WebGPU::initialize_bind_group_layout(wgpu::Device device) {
//     std::cout << "Preparing bind group layout..." << std::endl;
//     std::vector<wgpu::BindGroupLayoutEntry> bindings(3, wgpu::Default);

//     // atom buffers
//     bindings[0].binding = 0;
//     bindings[0].buffer.type = wgpu::BufferBindingType::ReadOnlyStorage;
//     bindings[0].visibility = wgpu::ShaderStage::Compute;

//     bindings[1].binding = 1;
//     bindings[1].buffer.type = wgpu::BufferBindingType::ReadOnlyStorage;
//     bindings[1].visibility = wgpu::ShaderStage::Compute;

//     // histogram buffer
//     bindings[2].binding = 2;
//     bindings[2].buffer.type = wgpu::BufferBindingType::Storage;
//     bindings[2].visibility = wgpu::ShaderStage::Compute;

//     wgpu::BindGroupLayoutDescriptor bindgroup_layout;
//     bindgroup_layout.entryCount = (uint32_t) bindings.size();
//     bindgroup_layout.entries = bindings.data();
//     return device.createBindGroupLayout(bindgroup_layout);
// }

// ComputePipelines WebGPU::initialize_compute_pipelines(wgpu::Device device, wgpu::BindGroupLayout bind_group_layout) {
//     wgpu::ShaderModule compute_shader_module = load_shader_module("executable/gpu_debug.wgsl", device);
//     std::cout << "Shader module loaded: " << compute_shader_module << std::endl;

//     auto ctr = [&] (std::string_view shader_name) {
//         std::cout << "Creating compute pipeline for shader: " << shader_name << std::endl;
//         wgpu::PipelineLayoutDescriptor pipeline_layout_desc;
//         pipeline_layout_desc.bindGroupLayoutCount = 1;
//         pipeline_layout_desc.bindGroupLayouts = (WGPUBindGroupLayout*) &bind_group_layout;
//         auto pipeline_layout = device.createPipelineLayout(pipeline_layout_desc);

//         wgpu::ComputePipelineDescriptor pipeline_descriptor = wgpu::Default;
//         pipeline_descriptor.compute.entryPoint.data = shader_name.data();
//         pipeline_descriptor.compute.entryPoint.length = shader_name.size();
//         pipeline_descriptor.compute.module = compute_shader_module;
//         pipeline_descriptor.layout = pipeline_layout;
//         return device.createComputePipeline(pipeline_descriptor);
//     };

//     return {
//         .self=ctr("calculate_self"),
//         .cross=ctr("calculate_cross")
//     };
// }

// BufferData create_buffers_test(wgpu::Device device) {
//     std::cout << "Creating buffers..." << std::endl;
//     static std::string atom_buffer_1_name = "atom buffer 1";
//     wgpu::BufferDescriptor atom_buffer_1_desc;
//     atom_buffer_1_desc.label.data = atom_buffer_1_name.data();
//     atom_buffer_1_desc.label.length = atom_buffer_1_name.size();
//     atom_buffer_1_desc.size = 8 * sizeof(float);
//     atom_buffer_1_desc.mappedAtCreation = false;
//     atom_buffer_1_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
//     wgpu::Buffer atom_buffer_1 = device.createBuffer(atom_buffer_1_desc);

//     static std::string atom_buffer_2_name = "atom buffer 2";
//     wgpu::BufferDescriptor atom_buffer_2_desc;
//     atom_buffer_2_desc.label.data = atom_buffer_2_name.data();
//     atom_buffer_2_desc.label.length = atom_buffer_2_name.size();
//     atom_buffer_2_desc.size = 8 * sizeof(float);
//     atom_buffer_2_desc.mappedAtCreation = false;
//     atom_buffer_2_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
//     wgpu::Buffer atom_buffer_2 = device.createBuffer(atom_buffer_2_desc);

//     static std::string histogram_buffer_name = "histogram buffer";
//     wgpu::BufferDescriptor histogram_buffer_desc;
//     histogram_buffer_desc.label.data = histogram_buffer_name.data();
//     histogram_buffer_desc.label.length = histogram_buffer_name.size();
//     histogram_buffer_desc.size = 8 * sizeof(float);
//     histogram_buffer_desc.mappedAtCreation = false;
//     histogram_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
//     wgpu::Buffer histogram_buffer = device.createBuffer(histogram_buffer_desc);
//     std::cout << "Buffers created." << std::endl;

//     return {atom_buffer_1, atom_buffer_2, histogram_buffer};
// }

// std::vector<float> convert_data(const std::vector<Atom>& atoms) {
//     std::vector<float> data(atoms.size()*4);
//     for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
//         const auto& atom = atoms[i];
//         data[i*4 + 0] = atom.x();
//         data[i*4 + 1] = atom.y();
//         data[i*4 + 2] = atom.z();
//         data[i*4 + 3] = atom.w;
//     }
//     return data;
// }

// wgpu::Buffer create_buffer(wgpu::Device device, const std::vector<Atom>& atoms) {
//     auto queue = device.getQueue();
//     auto data1 = convert_data(atoms);
//     wgpu::BufferDescriptor atom_buffer_desc;
//     atom_buffer_desc.size = data1.size()*sizeof(float);
//     atom_buffer_desc.mappedAtCreation = false;
//     atom_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
//     wgpu::Buffer atom_buffer = device.createBuffer(atom_buffer_desc);
//     queue.writeBuffer(atom_buffer, 0, data1.data(), atom_buffer.getSize());
//     return atom_buffer;
// }

// BufferData create_buffers(wgpu::Device device, const std::vector<Atom>& atoms) {
//     assert(false && "Not implemented yet");
//     auto atom_buffer_1 = create_buffer(device, atoms);

//     wgpu::BufferDescriptor histogram_buffer_desc;
//     histogram_buffer_desc.size = constants::axes::d_vals.size()*sizeof(float);
//     histogram_buffer_desc.mappedAtCreation = false;
//     histogram_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
//     wgpu::Buffer histogram_buffer = device.createBuffer(histogram_buffer_desc);
//     return {atom_buffer_1, {}, histogram_buffer};
// }

// BufferData create_buffers(wgpu::Device device, const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2) {
//     auto atom_buffer_1 = create_buffer(device, atoms1);
//     auto atom_buffer_2 = create_buffer(device, atoms2);

//     wgpu::BufferDescriptor histogram_buffer_desc;
//     histogram_buffer_desc.size = constants::axes::d_vals.size()*sizeof(float);
//     histogram_buffer_desc.mappedAtCreation = false;
//     histogram_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
//     wgpu::Buffer histogram_buffer = device.createBuffer(histogram_buffer_desc);
//     return {atom_buffer_1, atom_buffer_2, histogram_buffer};
// }

// wgpu::BindGroup WebGPU::assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer, wgpu::BindGroupLayout bind_group_layout) {
//     assert(buffer1 && "Buffer 1 is null");
//     assert(buffer2 && "Buffer 2 is null");
//     assert(histogram_buffer && "Histogram buffer is null");
//     std::cout << "Assigning buffers to bind group..." << std::endl;

//     std::vector<wgpu::BindGroupEntry> entries(3, wgpu::Default);
//     entries[0].binding = 0;
//     entries[0].buffer = buffer1;
//     entries[0].size = buffer1.getSize();

//     entries[1].binding = 1;
//     entries[1].buffer = buffer2;
//     entries[1].size = buffer2.getSize();

//     entries[2].binding = 2;
//     entries[2].buffer = histogram_buffer;
//     entries[2].size = histogram_buffer.getSize();

//     wgpu::BindGroupDescriptor bind_group_desc;
//     bind_group_desc.layout = bind_group_layout;
//     bind_group_desc.entryCount = entries.size();
//     bind_group_desc.entries = (WGPUBindGroupEntry*) entries.data();
//     return device.createBindGroup(bind_group_desc);
// }

// void WebGPU::test_fill_buffers(wgpu::Queue queue, wgpu::Buffer buffer1, wgpu::Buffer buffer2) {
//     std::cout << "Filling buffers..." << std::endl;
//     std::array<float, 8> data1 = {0, 0, 1, 2, 0, 0, 2, 4};
//     std::array<float, 8> data2 = {1, 0, 1, 6, 1, 0, 2, 8};
//     queue.writeBuffer(buffer1, 0, data1.data(), buffer1.getSize());
//     queue.writeBuffer(buffer2, 0, data2.data(), buffer2.getSize());
// }

// void WebGPU::submit_self(const std::vector<Atom>& atoms) {
//     std::cout << "Calculating self..." << std::endl;
//     wgpu::CommandEncoder encoder = device.createCommandEncoder();

//     buffers = create_buffers(device, atoms, atoms); //! use single buffer version for self calculation
//     wgpu::BindGroup bind_group = assign_buffers(device, buffers.atom_1, buffers.atom_2, buffers.histogram, bind_group_layout);

//     // submit work
//     wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
//     auto compute_pass = encoder.beginComputePass(compute_pass_desc);
//     compute_pass.setPipeline(pipelines.self);
//     compute_pass.setBindGroup(0, bind_group, 0, nullptr);
//     compute_pass.dispatchWorkgroups(1, 1, 1);
//     compute_pass.end();

//     // retrieve data
//     wgpu::BufferDescriptor readback_buffer_desc;
//     readback_buffer_desc.size = buffers.histogram.getSize();
//     readback_buffer_desc.mappedAtCreation = false;
//     readback_buffer_desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
//     wgpu::Buffer readback_buffer = device.createBuffer(readback_buffer_desc);

//     assert(readback_buffer.getSize() == buffers.histogram.getSize() && "Readback buffer size does not match histogram buffer size.");
//     encoder.copyBufferToBuffer(buffers.histogram, 0, readback_buffer, 0, readback_buffer.getSize());

//     wgpu::CommandBufferDescriptor command_buffer_desc;
//     command_buffer_desc.nextInChain = nullptr;
//     auto command_buffer = encoder.finish(command_buffer_desc);
//     encoder.release();

//     device.getQueue().submit(command_buffer);
//     command_buffer.release();
//     readback_buffers.push_back(readback_buffer);    
// }

// void WebGPU::submit_cross(const std::vector<Atom>& atoms1, const std::vector<Atom>& atoms2) {
//     std::cout << "Calculating cross..." << std::endl;
//     wgpu::CommandEncoder encoder = device.createCommandEncoder();

//     buffers = create_buffers(device, atoms1, atoms2);
//     wgpu::BindGroup bind_group = assign_buffers(device, buffers.atom_1, buffers.atom_2, buffers.histogram, bind_group_layout);

//     // submit work
//     wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
//     auto compute_pass = encoder.beginComputePass(compute_pass_desc);
//     compute_pass.setPipeline(pipelines.self);
//     compute_pass.setBindGroup(0, bind_group, 0, nullptr);
//     compute_pass.dispatchWorkgroups(1, 1, 1);
//     compute_pass.end();

//     // retrieve data
//     wgpu::BufferDescriptor readback_buffer_desc;
//     readback_buffer_desc.size = buffers.histogram.getSize();
//     readback_buffer_desc.mappedAtCreation = false;
//     readback_buffer_desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
//     wgpu::Buffer readback_buffer = device.createBuffer(readback_buffer_desc);
//     assert(readback_buffer.getSize() == buffers.histogram.getSize() && "Readback buffer size does not match histogram buffer size.");
//     encoder.copyBufferToBuffer(buffers.histogram, 0, readback_buffer, 0, readback_buffer.getSize());

//     wgpu::CommandBufferDescriptor command_buffer_desc;
//     command_buffer_desc.nextInChain = nullptr;
//     auto command_buffer = encoder.finish(command_buffer_desc);
//     encoder.release();

//     device.getQueue().submit(command_buffer);
//     command_buffer.release();
//     readback_buffers.push_back(readback_buffer);
// }

#include <gpu/WebGPU/WebGPUController.h>
int main(int, char const *[]) {
    std::vector<Atom> atoms = {
        Atom({-1, -1, -1}, 1), Atom({-1, 1, -1}, 1),
        Atom({ 1, -1, -1}, 1), Atom({ 1, 1, -1}, 1),
        Atom({-1, -1,  1}, 1), Atom({-1, 1,  1}, 1),
        Atom({ 1, -1,  1}, 1), Atom({ 1, 1,  1}, 1)
    };

    ausaxs::gpu::WebGPU webgpu;
    webgpu.submit_self(atoms);
    auto res = webgpu.run();

    // std::cout << "Results:" << std::endl;
    // for (int i = 0; i < 20; ++i) {
    //     std::cout << res[i] << std::endl;
    // }

    return 0;
}