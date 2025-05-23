#include <webgpu/webgpu.hpp>

#include <iostream>
#include <filesystem>
#include <fstream>

wgpu::ShaderModule load_shader_module(const std::filesystem::path& path, wgpu::Device device) {
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

wgpu::Adapter get_adapter() {
    // we do not need any special features, so we can use the default descriptor
    std::cout << "Creating WebGPU instance..." << std::endl;
    wgpu::InstanceDescriptor descriptor;
    wgpu::Instance instance = wgpu::createInstance(descriptor);

    if (!instance) {
        std::cout << "Could not create WebGPU instance." << std::endl;
        exit(0);
    } else {
        std::cout << "WebGPU instance created." << std::endl;
    }

    struct UserData {
        wgpu::Adapter adapter;
        bool request_finished = false;
    }; 
    UserData user_data;

    {   // perform the adapter request        
        wgpu::RequestAdapterOptions options;
        options.powerPreference = wgpu::PowerPreference::HighPerformance;

        // define callback function triggered when the adapter is ready (or failed)
        wgpu::RequestAdapterCallbackInfo adapter_callback;
        adapter_callback.userdata1 = &user_data;
        adapter_callback.mode = wgpu::CallbackMode::WaitAnyOnly;
        adapter_callback.callback = [] (WGPURequestAdapterStatus status, WGPUAdapter adapter, WGPUStringView message, void* user_data_p, void*) -> void {
            UserData& user_data = *reinterpret_cast<UserData*>(user_data_p);
            if (status == wgpu::RequestAdapterStatus::Success) {
                user_data.adapter = adapter;
            } else {
                std::cout << "Could not get WebGPU adapter: " << std::string(message.data, message.length) << std::endl;
            }
            user_data.request_finished = true;
        };
        instance.requestAdapter(options, adapter_callback);
    }

    assert(user_data.request_finished);

    // print adapter info
    wgpu::AdapterInfo info;
    user_data.adapter.getInfo(&info);
    std::cout << "WebGPU: Acquired device: " << std::string(info.device.data, info.device.length) << std::endl;

    return user_data.adapter;
}

wgpu::Device get_device() {
    // create a device with default descriptor
    auto adapter = get_adapter();
    wgpu::DeviceDescriptor device_descriptor;

    {   
        // setup callback function for disconnected devices
        wgpu::DeviceLostCallbackInfo device_callback;
        device_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
        device_callback.callback = [] (WGPUDeviceImpl *const *, WGPUDeviceLostReason reason, WGPUStringView message, void *, void *) -> void {
            std::cout << "Device lost: reason " << reason;
            if (message.data) std::cout << " (" << std::string(message.data, message.length) << ")";
            std::cout << std::endl;
        };

        // setup callback function for uncaptured errors
        wgpu::UncapturedErrorCallbackInfo error_callback;
        error_callback.callback = [] (WGPUDeviceImpl *const *, WGPUErrorType error, WGPUStringView message, void *, void *) -> void {
            std::cout << "Uncaptured error: " << error;
            if (message.data) std::cout << " (" << std::string(message.data, message.length) << ")";
            std::cout << std::endl;
        };

        device_descriptor.uncapturedErrorCallbackInfo = error_callback;
        device_descriptor.deviceLostCallbackInfo = device_callback;
    }

    struct UserData {
        WGPUDevice device = nullptr;
        bool requestEnded = false;
    };
    UserData userData;

    {   // perform the device request
        WGPURequestDeviceCallbackInfo device_callback;
        device_callback.userdata1 = &userData;
        device_callback.mode = wgpu::CallbackMode::WaitAnyOnly;
        device_callback.callback = [] (WGPURequestDeviceStatus status, WGPUDevice device, WGPUStringView message, void* user_data_p, void*) -> void {
            UserData& userData = *reinterpret_cast<UserData*>(user_data_p);
            if (status == wgpu::RequestDeviceStatus::Success) {
                userData.device = device;
            } else {
                std::cout << "Could not get WebGPU device: " << std::string(message.data, message.length) << std::endl;
            }
            userData.requestEnded = true;
        };
        adapter.requestDevice(device_descriptor, device_callback);    
    }

    assert(userData.requestEnded);
    adapter.release();
    return userData.device;
}

wgpu::BindGroupLayout initialize_bind_group_layout(wgpu::Device device) {
    std::cout << "Preparing bind group layout..." << std::endl;
    std::vector<wgpu::BindGroupLayoutEntry> bindings(3, wgpu::Default);

    // atom buffers
    bindings[0].binding = 0;
    bindings[0].buffer.type = wgpu::BufferBindingType::ReadOnlyStorage;
    bindings[0].visibility = wgpu::ShaderStage::Compute;
    bindings[1].binding = 1;
    bindings[1].buffer.type = wgpu::BufferBindingType::ReadOnlyStorage;
    bindings[1].visibility = wgpu::ShaderStage::Compute;

    // histogram buffer
    bindings[2].binding = 2;
    bindings[2].buffer.type = wgpu::BufferBindingType::Storage;
    bindings[2].visibility = wgpu::ShaderStage::Compute;

    wgpu::BindGroupLayoutDescriptor bindgroup_layout;
    bindgroup_layout.entryCount = (uint32_t) bindings.size();
    bindgroup_layout.entries = bindings.data();
    return device.createBindGroupLayout(bindgroup_layout);
}

wgpu::ComputePipeline initialize_compute_pipeline(wgpu::Device device, wgpu::BindGroupLayout bind_group_layout) {
    wgpu::ShaderModule compute_shader_module = load_shader_module("executable/gpu_debug.wgsl", device);
    std::cout << "Shader module loaded: " << compute_shader_module << std::endl;

    std::cout << "Creating compute pipeline descriptor..." << std::endl;
    std::string entry_point = "calculate_self";
    wgpu::ComputePipelineDescriptor pipeline_descriptor = wgpu::Default;
    pipeline_descriptor.compute.entryPoint.data = entry_point.data();
    pipeline_descriptor.compute.entryPoint.length = entry_point.size();
    pipeline_descriptor.compute.module = compute_shader_module;
    wgpu::ComputePipeline pipeline = device.createComputePipeline(pipeline_descriptor);

    std::cout << "Creating compute pipeline..." << std::endl;
    wgpu::PipelineLayoutDescriptor pipeline_layout_desc;
    pipeline_layout_desc.bindGroupLayoutCount = 1;
    pipeline_layout_desc.bindGroupLayouts = (WGPUBindGroupLayout*) &bind_group_layout;
    auto pipeline_layout = device.createPipelineLayout(pipeline_layout_desc);
    pipeline_descriptor.layout = pipeline_layout;
    return pipeline;
}

std::array<wgpu::Buffer, 3> create_buffers(wgpu::Device device) {
    std::cout << "Creating buffers..." << std::endl;
    wgpu::BufferDescriptor atom_buffer_1_desc;
    atom_buffer_1_desc.size = 8 * sizeof(float);
    atom_buffer_1_desc.mappedAtCreation = false;
    atom_buffer_1_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
    wgpu::Buffer atom_buffer_1 = device.createBuffer(atom_buffer_1_desc);

    wgpu::BufferDescriptor atom_buffer_2_desc;
    atom_buffer_2_desc.size = 8 * sizeof(float);
    atom_buffer_2_desc.mappedAtCreation = false;
    atom_buffer_2_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
    wgpu::Buffer atom_buffer_2 = device.createBuffer(atom_buffer_2_desc);

    wgpu::BufferDescriptor histogram_buffer_desc;
    histogram_buffer_desc.size = 8 * sizeof(float);
    histogram_buffer_desc.mappedAtCreation = false;
    histogram_buffer_desc.usage = wgpu::BufferUsage::Storage;
    wgpu::Buffer histogram_buffer = device.createBuffer(histogram_buffer_desc);
    std::cout << "Buffers created." << std::endl;

    return {atom_buffer_1, atom_buffer_2, histogram_buffer};
}

wgpu::BindGroup assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer, wgpu::BindGroupLayout bind_group_layout) {
    std::array<wgpu::BindGroupEntry, 3> entries;
    entries[0].binding = 0;
    entries[0].buffer = buffer1;
    entries[0].size = buffer1.getSize();
    entries[0].offset = 0;

    entries[1].binding = 1;
    entries[1].buffer = buffer2;
    entries[1].size = buffer2.getSize();
    entries[1].offset = 0;

    entries[2].binding = 2;
    entries[2].buffer = histogram_buffer;
    entries[2].size = histogram_buffer.getSize();
    entries[2].offset = 0;

    wgpu::BindGroupDescriptor bind_group_desc;
    bind_group_desc.layout = bind_group_layout;
    bind_group_desc.entryCount = entries.size();
    bind_group_desc.entries = entries.data();
    return device.createBindGroup(bind_group_desc);
}

int main(int, char const *[]) {
    auto device = get_device();
    auto bind_group_layout = initialize_bind_group_layout(device);
    auto pipeline = initialize_compute_pipeline(device, bind_group_layout);

    auto[atom_buffer_1, atom_buffer_2, histogram_buffer] = create_buffers(device);
    auto bind_group = assign_buffers(device, atom_buffer_1, atom_buffer_2, histogram_buffer, bind_group_layout);
    auto encoder = device.createCommandEncoder();

    // submit work
    wgpu::ComputePassDescriptor compute_pass_desc = wgpu::Default;
    auto compute_pass = encoder.beginComputePass(compute_pass_desc);
    compute_pass.setPipeline(pipeline);
    compute_pass.setBindGroup(0, bind_group, 0, nullptr);
    compute_pass.dispatchWorkgroups(1, 1, 1);
    compute_pass.end();

    // retrieve data
    wgpu::BufferDescriptor readback_buffer_desc;
    readback_buffer_desc.size = 8*sizeof(float);
    readback_buffer_desc.mappedAtCreation = false;
    readback_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
    wgpu::Buffer readback_buffer = device.createBuffer(readback_buffer_desc);
    encoder.copyBufferToBuffer(histogram_buffer, 0, readback_buffer, 0, readback_buffer_desc.size);

    // print
    bool done = false;
    wgpu::BufferMapCallbackInfo map_callback;
    map_callback.mode = wgpu::CallbackMode::WaitAnyOnly;
    map_callback.callback = [] (WGPUMapAsyncStatus, WGPUStringView, void* p_readback_buffer, void* p_done) -> void {
        wgpu::Buffer readback_buffer = *reinterpret_cast<wgpu::Buffer*>(p_readback_buffer);
        bool* done = reinterpret_cast<bool*>(p_done);
        const float* output = (const float*) readback_buffer.getConstMappedRange(0, readback_buffer.getSize());
        for (int i = 0; i < 8; ++i) {
            std::cout << "output " << i << ": " << output[i] << std::endl;
        }
        readback_buffer.unmap();
        *done = true;
    };

    while (!done) {
        #ifdef WEBGPU_BACKEND_WGPU
            device.getQueue().submit(0, nullptr);
        #else
            throw std::runtime_error("WebGPU backend not supported.");
        #endif
    }

    return 0;
}