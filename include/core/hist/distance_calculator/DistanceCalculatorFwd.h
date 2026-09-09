// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs::hist::distance_calculator {
    template<bool weighted_bins, bool variable_bin_width> class SimpleCalculator;
    template<bool weighted_bins, bool variable_bin_width> class SimpleCPU;
    template<bool weighted_bins, bool variable_bin_width> class SimpleGPU;
}
