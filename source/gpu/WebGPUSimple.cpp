#include <gpu/WebGPUSimple.h>

#define WEBGPU_CPP_IMPLEMENTATION
#include <webgpu/webgpu.hpp>

#include <filesystem>
#include <fstream>

using namespace ausaxs;
using namespace ausaxs::hist::distance_calculator;

template<bool weighted_bins>
WebGPUSimple<weighted_bins>::WebGPUSimple() {
    initialize();
}

namespace {
void print_adapter_info(wgpu::Adapter adapter) {
    wgpu::AdapterInfo info;
    adapter.getInfo(&info);
    std::cout << "WebGPU: Acquired device: " << std::string(info.device.data, info.device.length) << std::endl;
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
    print_adapter_info(user_data.adapter);
    return user_data.adapter;
}

wgpu::Device get_device(wgpu::Adapter adapter) {
    // create a device with default descriptor
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

wgpu::Queue get_queue(wgpu::Device device) {
    auto queue = device.getQueue();

    wgpu::QueueWorkDoneCallbackInfo queue_callback;
    queue_callback.mode = wgpu::CallbackMode::WaitAnyOnly;
    queue_callback.callback = [] (WGPUQueueWorkDoneStatus status, void *, void *) -> void {
        if (status == wgpu::QueueWorkDoneStatus::Success) {
            std::cout << "Queue finished successfully." << std::endl;
        } else {
            std::cout << "Queue failed its job." << std::endl;
        }
    };
    queue.onSubmittedWorkDone(queue_callback);
    return queue;
}

wgpu::ComputePassEncoder get_encoder(wgpu::Device device, wgpu::Queue queue) {
    wgpu::CommandEncoderDescriptor encoder_descriptor = wgpu::Default;
    wgpu::CommandEncoder encoder = device.createCommandEncoder(encoder_descriptor);

    // compute pass
    wgpu::ComputePassDescriptor compute_pass_descriptor = wgpu::Default;
    wgpu::ComputePassEncoder compute_pass = encoder.beginComputePass(compute_pass_descriptor);
    compute_pass.end();

    wgpu::CommandBufferDescriptor command_buffer_descriptor = wgpu::Default;
    wgpu::CommandBuffer commands = encoder.finish(command_buffer_descriptor);
    queue.submit(commands);

    // cleanup
    compute_pass.release();
    encoder.release();
    commands.release();
    queue.release();
    return compute_pass;
}

wgpu::ShaderModule load_shader_module(const std::filesystem::path& path, wgpu::Device device) {
    std::ifstream file(path);
    if (!file.is_open()) {
        return nullptr;
    }
    file.seekg(0, std::ios::end);
    size_t size = file.tellg();
    std::string source(size, ' ');
    file.seekg(0);
    file.read(source.data(), size);

    // wgpu::ShaderModuleGLSLDescriptor shader_code_descriptor{};
    // shader_code_descriptor.chain.next = nullptr;
    // shader_code_descriptor.chain.sType = wgpu::SType::ShaderSourceWGSL;
    // shader_code_descriptor.code.data = source.c_str();
    // shader_code_descriptor.code.length = source.size();

    wgpu::ShaderModuleDescriptor shaderDesc{};
    
    #ifdef WEBGPU_BACKEND_WGPU
        shader_code_descriptor.defineCount = 0;
        shader_code_descriptor.defines = nullptr;
    #endif
    
    // shaderDesc.nextInChain = &shader_code_descriptor.chain;
    return device.createShaderModule(shaderDesc);
}

wgpu::ComputePipeline prepare_compute_pipeline(wgpu::Device device) {
    std::cout << "Preparing compute pipeline..." << std::endl;
    wgpu::ShaderModule module;
    load_shader_module("source/gpu/WebGPUSimple.wgsl", device);

    std::string func = "calculate_self";
    wgpu::ComputePipelineDescriptor pipeline_descriptor = wgpu::Default;
    pipeline_descriptor.compute.module = module;
    pipeline_descriptor.compute.entryPoint.data = func.data();
    pipeline_descriptor.compute.entryPoint.length = func.size();
    return device.createComputePipeline(pipeline_descriptor);
}

wgpu::BindGroupLayout prepare_bindings(wgpu::Device device) {
    std::cout << "Preparing bindings..." << std::endl;
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

wgpu::Device device;
wgpu::Queue queue;
wgpu::BindGroupLayout bindings;
wgpu::ComputePipeline pipeline;
wgpu::Buffer atom_buffer_1;
wgpu::Buffer atom_buffer_2;
wgpu::Buffer histogram_buffer;

void fill_input_buffers(const hist::detail::CompactCoordinates& a) {
    wgpu::BufferDescriptor atom_buffer_desc;
    atom_buffer_desc.size = a.size() * sizeof(hist::detail::CompactCoordinatesData);
    atom_buffer_desc.mappedAtCreation = false;
    atom_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
    wgpu::Buffer atom_buffer_1 = device.createBuffer(atom_buffer_desc);

    wgpu::BufferDescriptor histogram_buffer_desc;
    histogram_buffer_desc.size = constants::axes::d_axis.bins * sizeof(constants::axes::d_type);
    histogram_buffer_desc.mappedAtCreation = false;
    histogram_buffer_desc.usage = wgpu::BufferUsage::Storage;
    wgpu::Buffer histogram_buffer = device.createBuffer(histogram_buffer_desc);

    queue.writeBuffer(atom_buffer_1, 0, a.get_data().data(), atom_buffer_desc.size);
}

void fill_input_buffers(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2) {
    wgpu::BufferDescriptor atom_buffer_1_desc;
    atom_buffer_1_desc.size = a1.size() * sizeof(hist::detail::CompactCoordinatesData);
    atom_buffer_1_desc.mappedAtCreation = false;
    atom_buffer_1_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
    wgpu::Buffer atom_buffer_1 = device.createBuffer(atom_buffer_1_desc);

    wgpu::BufferDescriptor atom_buffer_2_desc;
    atom_buffer_2_desc.size = a2.size() * sizeof(hist::detail::CompactCoordinatesData);
    atom_buffer_2_desc.mappedAtCreation = false;
    atom_buffer_2_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
    wgpu::Buffer atom_buffer_2 = device.createBuffer(atom_buffer_2_desc);

    wgpu::BufferDescriptor histogram_buffer_desc;
    histogram_buffer_desc.size = constants::axes::d_axis.bins * sizeof(constants::axes::d_type);
    histogram_buffer_desc.mappedAtCreation = false;
    histogram_buffer_desc.usage = wgpu::BufferUsage::Storage;
    wgpu::Buffer histogram_buffer = device.createBuffer(histogram_buffer_desc);

    queue.writeBuffer(atom_buffer_1, 0, a1.get_data().data(), atom_buffer_1_desc.size);
    queue.writeBuffer(atom_buffer_2, 0, a2.get_data().data(), atom_buffer_2_desc.size);
}
}

template<bool weighted_bins>
void WebGPUSimple<weighted_bins>::initialize() {
    device = get_device(get_adapter());
    queue = get_queue(device);
    bindings = prepare_bindings(device);
    pipeline = prepare_compute_pipeline(device);
}

template<bool weighted_bins>
WebGPUSimple<weighted_bins>::~WebGPUSimple() {
    device.release();
    queue.release();
}

template<bool weighted_bins>
int WebGPUSimple<weighted_bins>::enqueue_calculate_self(const hist::detail::CompactCoordinates& a, int, int) {
    fill_input_buffers(a);
    return 0;
}

template<bool weighted_bins>
int WebGPUSimple<weighted_bins>::enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2, int, int) {
    fill_input_buffers(a1, a2);
    return 0;
}

template<bool weighted_bins>
SimpleKernel<weighted_bins>::run_result WebGPUSimple<weighted_bins>::run() {
    return {};
}

template struct ausaxs::hist::distance_calculator::WebGPUSimple<false>;
template struct ausaxs::hist::distance_calculator::WebGPUSimple<true>;