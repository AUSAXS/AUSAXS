// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/WebGPU/GPUInstance.h>
#include <utility/Exceptions.h>
#include <utility/Logging.h>

#include <thread>

using namespace ausaxs;
using namespace ausaxs::gpu;

namespace {
    wgpu::Instance create_instance() {
        wgpu::InstanceDescriptor descriptor;
        wgpu::Instance instance = wgpu::createInstance(descriptor);
        if (!instance) {throw except::runtime_error("GPUInstance: could not create a WebGPU instance.");}
        return instance;
    }

    wgpu::Adapter get_adapter(wgpu::Instance instance) {
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
                    logging::log("GPUInstance: could not get a WebGPU adapter: " + std::string(message.data, message.length));
                }
                done = true;
            };
            instance.requestAdapter(options, adapter_callback);
            instance.processEvents();
        }

        if (!done || !adapter) {throw except::runtime_error("GPUInstance: could not get a WebGPU adapter.");}
        wgpu::AdapterInfo info;
        adapter.getInfo(&info);
        logging::log("GPUInstance: acquired device \"" + std::string(info.device.data, info.device.length) + "\".");
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
                std::string msg = "GPUInstance: device lost (reason " + std::to_string(static_cast<int>(reason)) + ")";
                if (message.data) {msg += ": " + std::string(message.data, message.length);}
                logging::log(msg);
            };

            // setup callback function for uncaptured errors
            wgpu::UncapturedErrorCallbackInfo error_callback;
            error_callback.callback = [] (WGPUDeviceImpl *const *, WGPUErrorType error, WGPUStringView message, void *, void *) -> void {
                std::string msg = "GPUInstance: uncaptured error (type " + std::to_string(static_cast<int>(error)) + ")";
                if (message.data) {msg += ": " + std::string(message.data, message.length);}
                logging::log(msg);
                std::cerr << msg << std::endl;
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
                    logging::log("GPUInstance: could not get a WebGPU device: " + std::string(message.data, message.length));
                }
                done = true;
            };
            adapter.requestDevice(device_descriptor, device_callback);
            instance.processEvents();
        }

        adapter.release();
        if (!done || !device) {throw except::runtime_error("GPUInstance: could not get a WebGPU device.");}
        return device;
    }
}

GPUInstance::GPUInstance() :
    instance(create_instance()),
    device(get_device(instance)),
    queue(device.getQueue())
{}

GPUInstance& GPUInstance::get() {
    static GPUInstance instance;
    return instance;
}

bool GPUInstance::available() {
    try {
        get();
        return true;
    } catch (const std::exception& e) {
        logging::log(std::string("GPUInstance::available: no usable device: ") + e.what());
        return false;
    }
}

void GPUInstance::process() {
    instance.processEvents();
}

void GPUInstance::wait(bool& done) {
    while (!done) {
        instance.processEvents();
        if (!done) {std::this_thread::yield();}
    }
}
