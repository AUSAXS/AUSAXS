// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::hist::distance_calculator {
    template<bool weighted_bins, bool variable_bin_width, bool unit_weights = false> class Calculator;
    template<bool weighted_bins> class HistogramStore;

    namespace detail {
        template<bool weighted_bins, bool variable_bin_width, bool unit_weights> class CalculatorCPU;
        template<bool weighted_bins, bool variable_bin_width, bool unit_weights> class GPUKernel;
    }
}
