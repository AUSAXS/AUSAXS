// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>

// Base template helpers - core distance evaluation functions used by all histogram managers
// These provide the low-level SIMD distance calculations

namespace ausaxs::detail::add8 {
    template<bool weighted_bins, bool variable_bin_width>
    inline auto evaluate(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        if constexpr (weighted_bins) {
            return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
        } else {
            return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
        }
    }
}

namespace ausaxs::detail::add4 {
    template<bool weighted_bins, bool variable_bin_width>
    inline auto evaluate(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        if constexpr (weighted_bins) {
            return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
        } else {
            return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
        }
    }
}

namespace ausaxs::detail::add1 {
    template<bool weighted_bins, bool variable_bin_width>
    inline auto evaluate(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        if constexpr (weighted_bins) {
            return data_i[i].evaluate(data_j[j]);
        } else {
            return data_i[i].evaluate_rounded(data_j[j]);
        }
    }
}

namespace ausaxs::detail::add8 {
    template<bool weighted_bins, bool variable_bin_width>
    inline auto evaluate(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_j, int i, int j) {
        if constexpr (weighted_bins) {
            return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
        } else {
            return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
        }
    }
}

namespace ausaxs::detail::add4 {
    template<bool weighted_bins, bool variable_bin_width>
    inline auto evaluate(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_j, int i, int j) {
        if constexpr (weighted_bins) {
            return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
        } else {
            return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
        }
    }
}

namespace ausaxs::detail::add1 {
    template<bool weighted_bins, bool variable_bin_width>
    inline auto evaluate(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_j, int i, int j) {
        if constexpr (weighted_bins) {
            return data_i[i].evaluate(data_j[j]);
        } else {
            return data_i[i].evaluate_rounded(data_j[j]);
        }
    }
}

// Mixed-type overloads for CompactCoordinatesFF + CompactCoordinates (water + excluded volume)
namespace ausaxs::detail::add8 {
    template<bool weighted_bins, bool variable_bin_width>
    inline auto evaluate(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        return evaluate<weighted_bins>(reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_width>&>(data_i), data_j, i, j);
    }
}

namespace ausaxs::detail::add4 {
    template<bool weighted_bins, bool variable_bin_width>
    inline auto evaluate(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        return evaluate<weighted_bins>(reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_width>&>(data_i), data_j, i, j);
    }
}

namespace ausaxs::detail::add1 {
    template<bool weighted_bins, bool variable_bin_width>
    inline auto evaluate(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        return evaluate<weighted_bins>(reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_width>&>(data_i), data_j, i, j);
    }
}
