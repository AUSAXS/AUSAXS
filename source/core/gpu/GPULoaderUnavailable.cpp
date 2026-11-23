// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/GPULoader.h>
#include <utility/Console.h>

using namespace ausaxs::gpu;

const GPULoader::Backend& GPULoader::get() {
    static const Backend backend = [] {
        ausaxs::console::print_warning("GPU support is unavailable in the WebAssembly build.");
        return Backend{};
    }();
    return backend;
}

bool GPULoader::available() {return static_cast<bool>(get());}

void GPULoader::report_failure(std::string_view action, abi::Status /*unused*/) {
    ausaxs::console::print_warning("GPU support is unavailable in the WebAssembly build: " + std::string(action));
}

std::string GPULoader::device_name() {return "none (GPU support is unavailable in the WebAssembly build)";}
std::string GPULoader::path() {return {};}
std::string_view GPULoader::library_name() {return "libausaxs_gpu_sycl.so";}

bool GPULoader::is_host_device(std::string_view id) {
    const std::string_view api = abi::host_api;
    return id.starts_with(api) && id.substr(api.size()).starts_with(':');
}

GPULoader::Probe GPULoader::probe(const std::string& /*unused*/, const DeviceChooser& /*unused*/) {
    Probe result;
    result.error = "GPU support is unavailable in the WebAssembly build";
    return result;
}
