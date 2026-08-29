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

#include <memory>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Factory for the pairwise distance histogram kernels.
     *
     * The GPU kernel is only available if the library was compiled with GPU support, and is only used
     * if it is explicitly requested through settings::general::gpu. Everything else uses the CPU kernel.
     */
    template<bool weighted_bins, bool variable_bin_width>
    struct SimpleCalculator {
        static std::unique_ptr<SimpleKernel<weighted_bins, variable_bin_width>> create() {
            if (settings::general::gpu) {
                #ifdef AUSAXS_GPU
                    return std::make_unique<WebGPUSimple<weighted_bins, variable_bin_width>>();
                #else
                    static bool warned = false;
                    if (!warned) {
                        warned = true;
                        console::print_warning(
                            "settings::general::gpu: this build has no GPU support; using the CPU. "
                            "Configure with -DGPU=ON to enable it."
                        );
                    }
                #endif
            }
            return std::make_unique<SimpleCPU<weighted_bins, variable_bin_width>>();
        }
    };
}
