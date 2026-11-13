// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/GPULoader.h>

#include <gpu/GPUInstaller.h>
#include <settings/GeneralSettings.h>
#include <settings/InternalState.h>
#include <utility/Console.h>
#include <utility/Logging.h>

#include <cassert>
#include <cstdint>
#include <cstdlib>
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
    std::string loaded_path;

    std::string see_troubleshooting() {
        return "\nSee " + std::string(GPULoader::troubleshooting_url) + " for common fixes.";
    }

    // set by report_failure() and never cleared, so a device that has failed once is not used again
    bool device_failed = false;

    /**
     * @brief Limit the AdaptiveCpp runtime to reporting errors, unless the user asked for something else.
     *
     * Its warnings are aimed at developers. Every user without an NVIDIA driver would e.g. be told that its CUDA plugin could not be loaded,
     * since it finds devices by loading the plugin for every backend, and the CUDA one needs the driver's own library. Must be called before
     * the backend is loaded, since the runtime reads it when it starts.
     */
    void quiet_runtime() {
        if (std::getenv("ACPP_DEBUG_LEVEL") != nullptr || std::getenv("HIPSYCL_DEBUG_LEVEL") != nullptr) {return;}
        #ifdef _WIN32
            _putenv_s("ACPP_DEBUG_LEVEL", "1");
        #else
            setenv("ACPP_DEBUG_LEVEL", "1", 0);
        #endif
    }

    #ifdef _WIN32
        constexpr const char* library_name = "ausaxs_gpu_sycl.dll";

        using library_handle = HMODULE;
        /**
         * @brief Open a library, looking for its dependencies beside it before anywhere else.
         *
         * A packaged backend carries its runtime in its own directory, which LoadLibraryA would only search if it were on PATH. This is the 
         * counterpart of the $ORIGIN RUNPATH a Linux package carries. A bare name, or a backend whose runtime only PATH can find, still gets 
         * the default search as a fallback.
         */
        library_handle open_library(const std::string& path) {
            quiet_runtime();
            if (std::filesystem::path file(path); file.has_parent_path()) {
                // LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR is only accepted together with an absolute path
                std::string absolute = std::filesystem::absolute(file).string();
                if (auto handle = LoadLibraryExA(absolute.c_str(), nullptr, LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS)) {
                    return handle;
                }
            }
            return LoadLibraryA(path.c_str());
        }
        void* find_symbol(library_handle handle, const char* name) {
            return reinterpret_cast<void*>(GetProcAddress(handle, name));
        }
        void close_library(library_handle handle) {FreeLibrary(handle);}
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

        using library_handle = void*;
        library_handle open_library(const std::string& path) {
            quiet_runtime();
            if (std::filesystem::path file(path); file.has_parent_path()) {
                std::error_code ec;
                if (auto real = std::filesystem::canonical(file, ec); !ec) {
                    return dlopen(real.c_str(), RTLD_NOW | RTLD_LOCAL);
                }
            }
            return dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
        }
        void* find_symbol(library_handle handle, const char* name) {return dlsym(handle, name);}
        void close_library(library_handle handle) {dlclose(handle);}
        std::string open_error() {
            const char* error = dlerror();
            return (error != nullptr) ? error : "unknown error";
        }

        std::filesystem::path own_directory() {
            Dl_info info;
            if (dladdr(reinterpret_cast<const void*>(&own_directory), &info) == 0 || (info.dli_fname == nullptr)) {return {};}
            return std::filesystem::path(info.dli_fname).parent_path();
        }
    #endif

    /**
     * @brief The paths to try, in order.
     *
     * An explicitly configured path is used alone: if the user named a backend, silently using a different one is worse than not using the 
     * GPU at all. Otherwise the backend is looked for beside the library, which is where a packaged installation puts it, then in the folder 
     * `ausaxs gpu install` downloads it to, and is finally left to the platform loader, which covers a development build and a system-wide 
     * installation.
     */
    std::vector<std::string> candidate_paths() {
        if (!settings::general::gpu_library->empty()) {return {settings::general::gpu_library};}

        std::vector<std::string> candidates;
        if (auto directory = own_directory(); !directory.empty()) {
            candidates.push_back((directory/library_name).string());
        }
        if (auto installed = installer::installed_library(); !installed.empty()) {
            candidates.push_back(installed.path());
        }
        candidates.emplace_back(library_name);
        return candidates;
    }

    /**
     * @brief Resolve every entry point, or none of them.
     */
    bool resolve(library_handle handle, GPULoader::Backend& backend) {
        constexpr const char* not_a_backend = "the backend is missing one or more entry points; it is probably not a GPU backend";
        auto* abi_version = reinterpret_cast<abi::abi_version_fn>(find_symbol(handle, abi::symbol_abi_version));
        if (abi_version == nullptr) {
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
        backend.device_count = reinterpret_cast<abi::device_count_fn>(find_symbol(handle, abi::symbol_device_count));
        backend.device_id = reinterpret_cast<abi::device_id_fn>(find_symbol(handle, abi::symbol_device_id));
        backend.select_device = reinterpret_cast<abi::select_device_fn>(find_symbol(handle, abi::symbol_select_device));
        backend.begin = reinterpret_cast<abi::begin_fn>(find_symbol(handle, abi::symbol_begin));
        backend.submit = reinterpret_cast<abi::submit_fn>(find_symbol(handle, abi::symbol_submit));
        backend.finish_unweighted = reinterpret_cast<abi::finish_unweighted_fn>(find_symbol(handle, abi::symbol_finish_unweighted));
        backend.finish_weighted = reinterpret_cast<abi::finish_weighted_fn>(find_symbol(handle, abi::symbol_finish_weighted));
        if (
            (backend.available != nullptr) && (backend.device_name != nullptr) && (backend.last_error != nullptr) && 
            (backend.device_count != nullptr) && (backend.device_id != nullptr) && (backend.select_device != nullptr) && (backend.begin != nullptr) && 
            (backend.submit != nullptr) && (backend.finish_unweighted != nullptr) && (backend.finish_weighted != nullptr)) 
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
            library_handle handle = open_library(path);
            if (handle == nullptr) {
                // the platform error already names the file it failed to open
                attempts += "\n    " + open_error();
                continue;
            }

            logging::log("GPULoader: found GPU backend at " + path + ", attempting load");
            if (resolve(handle, backend)) {
                loaded_path = path;
                return backend;
            }
            logging::log("GPULoader: " + load_error);
            return {}; // note: intentionally leaks the handle
        }

        load_error = "no GPU backend could be opened:" + attempts;
        logging::log("GPULoader: " + load_error);
        return backend;
    }

    /**
     * @brief Identifiers of every GPU @a backend found, in its ranked order.
     *
     * The backend's host CPU device is left out unless settings::internal_state::gpu_expose_cpu is set. It only exists to test the backend on
     * machines without a GPU, and would be much slower than the CPU code for anyone who picked it by mistake.
     */
    std::vector<std::string> list_devices(const GPULoader::Backend& backend) {
        std::vector<std::string> devices;
        for (std::int32_t i = 0; i < backend.device_count(); ++i) {
            const char* id = backend.device_id(i);
            if (id == nullptr) {continue;}
            if (!settings::internal_state::gpu_expose_cpu && GPULoader::is_host_device(id)) {continue;}
            devices.emplace_back(id);
        }
        return devices;
    }

    /**
     * @brief Point @a backend at the device chosen when it was installed, if it is the installed one and a device was chosen.
     *
     * @return False if the chosen device could not be selected, in which case a warning has been printed.
     */
    bool select_installed_device(const GPULoader::Backend& backend) {
        auto installed = installer::installed_library();
        if (installed.empty() || loaded_path != installed.path()) {return true;}
        auto device = installer::installed_device();
        if (device.empty()) {return true;}

        if (backend.select_device(device.c_str()) == abi::Status::ok) {return true;}
        std::string found;
        for (const auto& id : list_devices(backend)) {found += "\n    " + id;}
        console::print_warning(
            "the GPU chosen when the backend was installed (" + device + ") is not present, so the GPU is not used. "
            "Run \"ausaxs gpu install\" again to choose another. The GPUs present are:" + (found.empty() ? std::string(" none") : found) + 
            see_troubleshooting()
        );
        return false;
    }
}

const GPULoader::Backend& GPULoader::get() {
    static const Backend backend = [] {
        Backend result = open();
        if (!result) {console::print_warning("no usable GPU backend: " + load_error + see_troubleshooting());}
        else if (!select_installed_device(result)) {
            device_failed = true; // the backend would otherwise go ahead on its own first choice
        } else {
            // resolve() only keeps a Backend whose every entry point resolved, and zeroes it otherwise, so a truthy one has them all
            assert(result.available != nullptr && result.device_name != nullptr && "GPULoader: an opened backend resolves every entry point.");
            if (!result.available()) {console::print_warning("no usable GPU backend: no usable device" + see_troubleshooting());}
            else {console::print_info("Using GPU backend: " + std::string(result.device_name()));}
        }
        return result;
    }();
    return backend;
}

bool GPULoader::available() {
    const auto& backend = get(); // first, since opening the backend may itself take the device out of use
    if (device_failed) {return false;}
    return backend && backend.available();
}

void GPULoader::report_failure(std::string_view action, abi::Status status) {
    if (device_failed) {return;}
    device_failed = true;

    const auto& backend = get();
    std::string reason = backend ? backend.last_error() : "no backend loaded";
    console::print_warning(
        "the GPU backend failed while trying to " + std::string(action) + " (status " + std::to_string(static_cast<int>(status)) + "): " 
        + reason + ". The GPU is not used again in this run." + see_troubleshooting()
    );
}

GPULoader::Probe GPULoader::probe(const std::string& path, const DeviceChooser& choose) {
    Probe result;
    library_handle handle = open_library(path);
    if (handle == nullptr) {
        result.error = open_error();
        return result;
    }

    Backend backend;
    if (!resolve(handle, backend)) {
        close_library(handle);
        result.error = load_error;
        return result;
    }
    result.loaded = true;

    result.devices = list_devices(backend);
    if (!result.devices.empty()) {
        std::size_t index = 0;
        if (auto choice = choose ? choose(result.devices) : std::nullopt) {
            result.chosen = true;
            if (*choice == no_device) {return result;}
            assert(*choice < result.devices.size() && "GPULoader::probe: the chooser must return the index of one of the devices.");
            index = *choice;
        }
        result.selected = result.devices[index];
        if (backend.select_device(result.selected.c_str()) != abi::Status::ok) {
            result.error = backend.last_error();
            return result;
        }
    }

    result.available = backend.available();
    if (result.available) {result.device = backend.device_name();}
    else {result.error = backend.last_error();}
    return result;
}

std::string GPULoader::path() {
    get();
    return loaded_path;
}

std::string_view GPULoader::library_name() {
    return ::library_name;
}

bool GPULoader::is_host_device(std::string_view id) {
    const std::string_view api = abi::host_api;
    return id.starts_with(api) && id.substr(api.size()).starts_with(':');
}

std::string GPULoader::device_name() {
    const auto& backend = get();
    if (!backend) {return "none (no backend loaded)";}
    if (device_failed) {return "none (the device is not in use)";}
    assert(backend.device_name != nullptr && "GPULoader: an opened backend resolves every entry point.");
    return backend.device_name();
}