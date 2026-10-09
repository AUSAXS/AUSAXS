// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/GPUBackendABI.h>

#include <string>
#include <string_view>

namespace ausaxs::gpu {
    /**
     * @brief Runtime access to a GPU backend. Instantiating anything in this namespace may trigger a just-in-time compilation of the kernels. 
     */
    class GPULoader {
        public:
            // @brief A backend, if one is loaded. 
            struct Backend {
                abi::available_fn available = nullptr;
                abi::device_name_fn device_name = nullptr;
                abi::last_error_fn last_error = nullptr;
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
                bool loaded = false;    // the library opened and implements this library's ABI
                bool available = false; // and it reports a usable device
                std::string device;     // the device's name, if available
                std::string error;      // why not, otherwise
            };

            /**
             * @brief Open the backend at @a path and ask it for a device, without making it the backend in use.
             *
             * For checking an installation. Note that a broken AdaptiveCpp runtime may terminate the process from inside this call, which cannot 
             * be caught; a caller that records an installation as good should only do so once this has returned. 
             */
            static Probe probe(const std::string& path);
    };
}
