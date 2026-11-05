// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/GPUBackendABI.h>

#include <cstddef>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace ausaxs::gpu {
    /**
     * @brief Runtime access to a GPU backend. Instantiating anything in this namespace may trigger a just-in-time compilation of the kernels. 
     */
    class GPULoader {
        public:
            // @brief Where the common fixes for a GPU that is not used are listed, for every message that reports one.
            static constexpr std::string_view troubleshooting_url = "https://github.com/AUSAXS/AUSAXS/wiki/Troubleshooting";

            // @brief A backend, if one is loaded.
            struct Backend {
                abi::available_fn available = nullptr;
                abi::device_name_fn device_name = nullptr;
                abi::last_error_fn last_error = nullptr;
                abi::device_count_fn device_count = nullptr;
                abi::device_id_fn device_id = nullptr;
                abi::select_device_fn select_device = nullptr;
                abi::begin_fn begin = nullptr;
                abi::submit_fn submit = nullptr;
                abi::finish_unweighted_fn finish_unweighted = nullptr;
                abi::finish_weighted_fn finish_weighted = nullptr;

                explicit operator bool() const {return begin != nullptr;}
            };

            /**
             * @brief Get the backend, opening it if this is the first call.
             *
             * The first call reports the outcome once: either the name of the device now in use, or a warning describing why there is none. Callers 
             * that only need to know whether to fall back to the CPU can just use available() and do not need to report anything themselves.
             * 
             * If the backend is the one `ausaxs gpu install` installed, the device chosen during the installation is used. If that device is no
             * longer present, the GPU is not used at all, rather than silently replaced by another one, which may well be slower than the CPU.
             */
            static const Backend& get();

            /**
             * @brief Whether a backend is loaded and reports a usable device.
             *        False for the rest of the process once report_failure() has been called.
             */
            static bool available();

            /**
             * @brief Report a backend call that failed, and take the device out of use.
             *
             * Prints one warning naming the first failure, and makes available() false from then on, so every later caller goes straight to the CPU. 
             *
             * @param action What the caller was trying to do, named in the warning. Format is "failed while trying to <action>".
             * @param status The status the backend returned.
             */
            static void report_failure(std::string_view action, abi::Status status);

            /**
             * @brief Name of the device in use, or a short description of why there is none.
             */
            static std::string device_name();

            /**
             * @brief Path of the loaded backend, or empty if none is loaded.
             */
            static std::string path();

            /**
             * @brief File name of the backend library on this platform.
             */
            static std::string_view library_name();

            // @brief The outcome of probe().
            struct Probe {
                bool loaded = false;                // the library opened and implements this library's ABI
                bool available = false;             // and it reports a usable device
                std::vector<std::string> devices;   // identifiers of every GPU the backend found, in its ranked order
                std::string selected;               // the identifier of the one that was tested, if any
                std::string device;                 // the device's name, if available
                std::string error;                  // why not, otherwise
            };

            /**
             * @brief Picks one of several GPUs, given their identifiers in the backend's ranked order. Returns the index of the chosen one.
             */
            using DeviceChooser = std::function<std::size_t(const std::vector<std::string>& devices)>;

            /**
             * @brief Open the backend at @a path and ask it for a device, without making it the backend in use.
             *
             * @param choose Asked to pick the device to test when the backend finds more than one. Without it, the first one is tested.
             *
             * For checking an installation. Note that a broken AdaptiveCpp runtime may terminate the process from inside this call, which cannot 
             * be caught; a caller that records an installation as good should only do so once this has returned. 
             */
            static Probe probe(const std::string& path, const DeviceChooser& choose = {});
    };
}
