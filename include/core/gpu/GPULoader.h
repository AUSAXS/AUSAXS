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
    };
}
