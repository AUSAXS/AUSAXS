// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::hist::distance_calculator {
    template<bool weighted_bins, bool variable_bin_width> class Calculator;
    template<bool weighted_bins, bool variable_bin_width> class CalculatorFF;
    template<bool weighted_bins> class HistogramStore;

    namespace detail {
        template<bool weighted_bins, bool variable_bin_width> class CalculatorCPU;
        template<bool weighted_bins, bool variable_bin_width> class GPUKernel;
    }
}
