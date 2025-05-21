#include <gpu/WebGPUSimple.h>

#define WEBGPU_CPP_IMPLEMENTATION
#include <webgpu/webgpu.hpp>

using namespace ausaxs;
using namespace ausaxs::hist::distance_calculator;

template<bool weighted_bins>
WebGPUSimple<weighted_bins>::WebGPUSimple() {
    initialize();
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
            std::cout << "Adapter callback triggered." << std::endl;
            UserData& user_data = *reinterpret_cast<UserData*>(user_data_p);
            if (status == wgpu::RequestAdapterStatus::Success) {
                std::cout << "Got WebGPU adapter: " << adapter << std::endl;
                user_data.adapter = adapter;
            } else {
                std::cout << "Could not get WebGPU adapter: " << std::string(message.data, message.length) << std::endl;
            }
            user_data.request_finished = true;
        };
        instance.requestAdapter(options, adapter_callback);
    }

    assert(user_data.request_finished);
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
            std::cout << "Device callback triggered." << std::endl;
            UserData& userData = *reinterpret_cast<UserData*>(user_data_p);
            if (status == wgpu::RequestDeviceStatus::Success) {
                std::cout << "Got WebGPU device: " << device << std::endl;
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

template<bool weighted_bins>
void WebGPUSimple<weighted_bins>::initialize() {
    auto adapter = get_adapter();
    auto device = get_device(adapter);
}

template<bool weighted_bins>
int WebGPUSimple<weighted_bins>::enqueue_calculate_self(const hist::detail::CompactCoordinates&, int, int) {
    return 0;
}

template<bool weighted_bins>
int WebGPUSimple<weighted_bins>::enqueue_calculate_cross(const hist::detail::CompactCoordinates&, const hist::detail::CompactCoordinates&, int, int) {
    return 0;
}

template<bool weighted_bins>
SimpleKernel<weighted_bins>::run_result WebGPUSimple<weighted_bins>::run() {
    return {};
}

template struct ausaxs::hist::distance_calculator::WebGPUSimple<false>;
template struct ausaxs::hist::distance_calculator::WebGPUSimple<true>;