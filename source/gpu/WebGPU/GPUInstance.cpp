#include <gpu/WebGPU/GPUInstance.h>

using namespace ausaxs;
using namespace ausaxs::gpu;

wgpu::Instance create_instance() {
    wgpu::InstanceDescriptor descriptor;
    wgpu::Instance instance = wgpu::createInstance(descriptor);
    assert(instance && "Could not create WebGPU instance.");
    return instance;
}

wgpu::Adapter get_adapter(wgpu::Instance instance) {
    if (!instance) {
        std::cout << "Could not create WebGPU instance." << std::endl;
        exit(0);
    } else {
        std::cout << "WebGPU instance created." << std::endl;
    }

    wgpu::Adapter adapter;
    bool done = false;
    {   // perform the adapter request
        wgpu::RequestAdapterOptions options;
        options.powerPreference = wgpu::PowerPreference::HighPerformance;

        wgpu::RequestAdapterCallbackInfo adapter_callback;
        adapter_callback.userdata1 = &adapter;
        adapter_callback.userdata2 = &done;
        adapter_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
        adapter_callback.callback = [] (WGPURequestAdapterStatus status, WGPUAdapter acquired_adapter, WGPUStringView message, void* adapter_p, void* done_p) -> void {
            wgpu::Adapter& adapter = *reinterpret_cast<wgpu::Adapter*>(adapter_p);
            bool& done = *reinterpret_cast<bool*>(done_p);
            if (status == wgpu::RequestAdapterStatus::Success) {
                adapter = acquired_adapter;
            } else {
                throw std::runtime_error("Could not get WebGPU adapter: " + std::string(message.data, message.length));
            }
            done = true;
        };
        instance.requestAdapter(options, adapter_callback);
        instance.processEvents();
    }

    assert(done);
    wgpu::AdapterInfo info;
    adapter.getInfo(&info);
    std::cout << "WebGPU: Acquired device: " << std::string(info.device.data, info.device.length) << std::endl;
    return adapter;
}

wgpu::Device get_device(wgpu::Instance instance) {
    auto adapter = get_adapter(instance);
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

    wgpu::Device device;
    bool done = false;
    {   // perform the device request
        WGPURequestDeviceCallbackInfo device_callback;
        device_callback.userdata1 = &device;
        device_callback.userdata2 = &done;
        device_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
        device_callback.callback = [] (WGPURequestDeviceStatus status, WGPUDevice acquired_device, WGPUStringView message, void* device_p, void* done_p) -> void {
            wgpu::Device& device = *reinterpret_cast<wgpu::Device*>(device_p);
            bool& done = *reinterpret_cast<bool*>(done_p);
            if (status == wgpu::RequestDeviceStatus::Success) {
                device = acquired_device;
            } else {
                throw std::runtime_error("Could not get WebGPU device: " + std::string(message.data, message.length));
            }
            done = true;
        };
        adapter.requestDevice(device_descriptor, device_callback);
        instance.processEvents();
    }

    assert(done);
    adapter.release();
    return device;
}

GPUInstance::GPUInstance() : 
    instance(create_instance()), 
    device(get_device(instance))
{}

void GPUInstance::process() {
    instance.processEvents();
}