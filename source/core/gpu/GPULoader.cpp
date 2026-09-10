// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/GPULoader.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Logging.h>

#include <filesystem>
#include <vector>

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#else
    #include <dlfcn.h>
#endif

using namespace ausaxs;
using namespace ausaxs::gpu;

namespace {
    std::string load_error;

    // set by report_failure() and never cleared, so a device that has failed once is not used again
    bool device_failed = false;

    #ifdef _WIN32
        constexpr const char* library_name = "ausaxs_gpu_sycl.dll";

        using handle_t = HMODULE;
        handle_t open_library(const std::string& path) {return LoadLibraryA(path.c_str());}
        void* find_symbol(handle_t handle, const char* name) {
            return reinterpret_cast<void*>(GetProcAddress(handle, name));
        }
        std::string open_error() {return "error " + std::to_string(GetLastError());}

        /// @brief Directory of the binary this code is part of.
        std::filesystem::path own_directory() {
            HMODULE self = nullptr;
            if (!GetModuleHandleExA(
                GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                reinterpret_cast<LPCSTR>(&own_directory), &self
            )) {return {};}
            char buffer[MAX_PATH];
            DWORD length = GetModuleFileNameA(self, buffer, MAX_PATH);
            if (length == 0 || length == MAX_PATH) {return {};}
            return std::filesystem::path(std::string(buffer, length)).parent_path();
        }
    #else
        #ifdef __APPLE__
            constexpr const char* library_name = "libausaxs_gpu_sycl.dylib";
        #else
            constexpr const char* library_name = "libausaxs_gpu_sycl.so";
        #endif

        using handle_t = void*;
        handle_t open_library(const std::string& path) {return dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);}
        void* find_symbol(handle_t handle, const char* name) {return dlsym(handle, name);}
        std::string open_error() {
            const char* error = dlerror();
            return error ? error : "unknown error";
        }

        std::filesystem::path own_directory() {
            Dl_info info;
            if (dladdr(reinterpret_cast<const void*>(&own_directory), &info) == 0 || !info.dli_fname) {return {};}
            return std::filesystem::path(info.dli_fname).parent_path();
        }
    #endif

    /**
     * @brief The paths to try, in order.
     *
     * An explicitly configured path is used alone: if the user named a backend, silently using a different one is worse than not using the 
     * GPU at all. Otherwise the backend is looked for beside the library, which is where a packaged installation puts it, and then left to 
     * the platform loader, which covers a development build and a system-wide installation.
     */
    std::vector<std::string> candidate_paths() {
        if (!settings::general::gpu_library.empty()) {return {settings::general::gpu_library};}

        std::vector<std::string> candidates;
        if (auto directory = own_directory(); !directory.empty()) {
            candidates.push_back((directory/library_name).string());
        }
        candidates.emplace_back(library_name);
        return candidates;
    }

    /**
     * @brief Resolve every entry point, or none of them.
     */
    bool resolve(handle_t handle, GPULoader::Backend& backend) {
        constexpr const char* not_a_backend = "the backend is missing one or more entry points; it is probably not a GPU backend";
        auto* abi_version = reinterpret_cast<abi::abi_version_fn>(find_symbol(handle, abi::symbol_abi_version));
        if (!abi_version) {
            load_error = not_a_backend;
            return false;
        }
        if (auto found = abi_version(); found != abi::version) {
            load_error =
                "the backend implements ABI version " + std::to_string(found) + ", but this library "
                "expects version " + std::to_string(abi::version) + "; they are from different releases";
            return false;
        }

        backend.available = reinterpret_cast<abi::available_fn>(find_symbol(handle, abi::symbol_available));
        backend.device_name = reinterpret_cast<abi::device_name_fn>(find_symbol(handle, abi::symbol_device_name));
        backend.last_error = reinterpret_cast<abi::last_error_fn>(find_symbol(handle, abi::symbol_last_error));
        backend.begin = reinterpret_cast<abi::begin_fn>(find_symbol(handle, abi::symbol_begin));
        backend.submit = reinterpret_cast<abi::submit_fn>(find_symbol(handle, abi::symbol_submit));
        backend.finish_unweighted = reinterpret_cast<abi::finish_unweighted_fn>(find_symbol(handle, abi::symbol_finish_unweighted));
        backend.finish_weighted = reinterpret_cast<abi::finish_weighted_fn>(find_symbol(handle, abi::symbol_finish_weighted));
        if (
            backend.available && backend.device_name && backend.last_error && backend.begin && 
            backend.submit && backend.finish_unweighted && backend.finish_weighted) 
        {
            return true;
        }

        load_error = not_a_backend;
        backend = {};
        return false;
    }

    GPULoader::Backend open() {
        GPULoader::Backend backend;
        std::string attempts;
        for (const auto& path : candidate_paths()) {
            handle_t handle = open_library(path);
            if (!handle) {
                // the platform error already names the file it failed to open
                attempts += "\n    " + open_error();
                continue;
            }

            logging::log("GPULoader: found GPU backend at " + path + ", attempting load");
            if (resolve(handle, backend)) {return backend;}
            logging::log("GPULoader: " + load_error);
            return {}; // note: intentionally leaks the handle
        }

        load_error = "no GPU backend could be opened:" + attempts;
        logging::log("GPULoader: " + load_error);
        return backend;
    }
}

const GPULoader::Backend& GPULoader::get() {
    static const Backend backend = [] {
        Backend result = open();
        if (!result) {console::print_warning("no usable GPU backend: " + load_error);}
        else if (!result.available()) {console::print_warning("no usable GPU backend: no usable device");}
        else {console::print_info("Using GPU backend: " + std::string(result.device_name()));}
        return result;
    }();
    return backend;
}

bool GPULoader::available() {
    if (device_failed) {return false;}
    const auto& backend = get();
    return backend && backend.available();
}

void GPULoader::report_failure(std::string_view action, abi::Status status) {
    if (device_failed) {return;}
    device_failed = true;

    const auto& backend = get();
    std::string reason = backend ? backend.last_error() : "no backend loaded";
    console::print_warning(
        "the GPU backend failed while trying to " + std::string(action) + " (status " + std::to_string(static_cast<int>(status)) + "): " 
        + reason + ". The GPU is not used again in this run."
    );
}

std::string GPULoader::device_name() {
    const auto& backend = get();
    if (!backend) {return "none (no backend loaded)";}
    return backend.device_name();
}