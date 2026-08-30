// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/SimpleKernel.h>
#include <hist/distance_calculator/SimpleCPU.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>

#ifdef AUSAXS_GPU
    #include <gpu/WebGPUSimple.h>
#endif

#ifdef AUSAXS_GPU_SYCL
    #include <gpu/SYCLSimple.h>
#endif

#include <memory>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Factory for the pairwise distance histogram kernels.
     *
     * A GPU kernel is only available if the library was compiled with support for one, and is only used
     * if it is explicitly requested through settings::general::gpu. Everything else uses the CPU kernel.
     * Which of the GPU kernels to use is chosen by settings::general::gpu_backend; if that one was not
     * compiled in, the other is used instead. Both defer the choice between the CPU and the GPU until
     * they know how large the calculation is, so this only picks what may be used, not what will be.
     */
    template<bool weighted_bins, bool variable_bin_width>
    struct SimpleCalculator {
        static std::unique_ptr<SimpleKernel<weighted_bins, variable_bin_width>> create() {
            if (settings::general::gpu) {
                #if defined(AUSAXS_GPU_SYCL) && defined(AUSAXS_GPU)
                    if (settings::general::gpu_backend == settings::general::GPUBackend::WEBGPU) {
                        return std::make_unique<WebGPUSimple<weighted_bins, variable_bin_width>>();
                    }
                    return std::make_unique<SYCLSimple<weighted_bins, variable_bin_width>>();
                #elif defined(AUSAXS_GPU_SYCL)
                    warn_once(settings::general::gpu_backend != settings::general::GPUBackend::SYCL,
                        "settings::general::gpu_backend: this build only has SYCL support; using that instead. "
                        "Configure with -DGPU=ON to enable the WebGPU backend."
                    );
                    return std::make_unique<SYCLSimple<weighted_bins, variable_bin_width>>();
                #elif defined(AUSAXS_GPU)
                    warn_once(settings::general::gpu_backend != settings::general::GPUBackend::WEBGPU,
                        "settings::general::gpu_backend: this build only has WebGPU support; using that instead. "
                        "Configure with -DSYCL=ON to enable the SYCL backend."
                    );
                    return std::make_unique<WebGPUSimple<weighted_bins, variable_bin_width>>();
                #else
                    warn_once(true,
                        "settings::general::gpu: this build has no GPU support; using the CPU. "
                        "Configure with -DSYCL=ON or -DGPU=ON to enable it."
                    );
                #endif
            }
            return std::make_unique<SimpleCPU<weighted_bins, variable_bin_width>>();
        }

        private:
            static void warn_once([[maybe_unused]] bool condition, [[maybe_unused]] std::string_view message) {
                static bool warned = false;
                if (condition && !warned) {
                    warned = true;
                    console::print_warning(std::string(message));
                }
            }
    };
}
