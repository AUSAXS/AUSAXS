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

    inline wgpu::Instance WebGPU::create_instance() {
        wgpu::InstanceDescriptor descriptor;
        wgpu::Instance instance = wgpu::createInstance(descriptor);
        assert(instance && "Could not create WebGPU instance.");
        return instance;
    }

    inline wgpu::Adapter WebGPU::get_adapter(wgpu::Instance instance) {
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
            adapter_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
            adapter_callback.callback = [] (WGPURequestAdapterStatus status, WGPUAdapter adapter, WGPUStringView message, void* user_data_p, void*) -> void {
                UserData& user_data = *reinterpret_cast<UserData*>(user_data_p);
                if (status == wgpu::RequestAdapterStatus::Success) {
                    std::cout << "Adapter found." << std::endl;
                    user_data.adapter = adapter;
                } else {
                    std::cout << "Could not get WebGPU adapter: " << std::string(message.data, message.length) << std::endl;
                }
                user_data.request_finished = true;
            };
            instance.requestAdapter(options, adapter_callback);
            instance.processEvents();
        }

        assert(user_data.request_finished);

        // print adapter info
        wgpu::AdapterInfo info;
        user_data.adapter.getInfo(&info);
        std::cout << "WebGPU: Acquired device: " << std::string(info.device.data, info.device.length) << std::endl;

        return user_data.adapter;
    }

    inline wgpu::Device WebGPU::get_device(wgpu::Instance instance) {
        // create a device with default descriptor
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

        struct UserData {
            WGPUDevice device = nullptr;
            bool requestEnded = false;
        };
        UserData userData;

        {   // perform the device request
            WGPURequestDeviceCallbackInfo device_callback;
            device_callback.userdata1 = &userData;
            device_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
            device_callback.callback = [] (WGPURequestDeviceStatus status, WGPUDevice device, WGPUStringView message, void* user_data_p, void*) -> void {
                UserData& userData = *reinterpret_cast<UserData*>(user_data_p);
                if (status == wgpu::RequestDeviceStatus::Success) {
                    std::cout << "Device found." << std::endl;
                    userData.device = device;
                } else {
                    std::cout << "Could not get WebGPU device: " << std::string(message.data, message.length) << std::endl;
                }
                userData.requestEnded = true;
            };
            adapter.requestDevice(device_descriptor, device_callback);
            instance.processEvents();
        }

        assert(userData.requestEnded);
        adapter.release();
        return userData.device;
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