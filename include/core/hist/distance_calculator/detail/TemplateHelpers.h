// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <form_factor/FormFactorType.h>

namespace ausaxs {
    constexpr int exv_bin = static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME);
}

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
        // Cast to base and delegate
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

namespace ausaxs {
    /**
     * @brief Calculate the distances between eight atoms and add them to the histogram. No excluded volume bin is added.
     * 
     * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
     * @tparam factor A multiplicative factor for the atomic weights. 
     * @param p The histogram to add the distances to.
     * @param data_i The first atom.
     * @param data_j The second atom.
     * @param i The index of the first atom.
     * @param j The index of the second atom.
     */
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinates<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add8::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 8; ++k) {
            p.template add<factor>(res.distances[k], res.weights[k]);
        }
    }

    /**
     * @brief Calculate the distances between four atoms and add them to the histogram. No excluded volume bin is added.
     * 
     * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
     * @tparam factor A multiplicative factor for the atomic weights. 
     * @param p The histogram to add the distances to.
     * @param data_i The first atom.
     * @param data_j The second atom.
     * @param i The index of the first atom.
     * @param j The index of the second atom.
     */
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinates<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add4::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 4; ++k) {
            p.template add<factor>(res.distances[k], res.weights[k]);
        }
    }

    /**
     * @brief Calculate the distances between two atoms and add them to the histogram. No excluded volume bin is added.
     * 
     * @tparam use_weighted_distribution Whether to keep track of the distances added to the bins. This is useful for weighting the bins later.
     * @tparam factor A multiplicative factor for the atomic weights. 
     * @param p The histogram to add the distances to.
     * @param data_i The first atom.
     * @param data_j The second atom.
     * @param i The index of the first atom.
     * @param j The index of the second atom.
     */
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinates<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add1::evaluate<weighted_bins>(data_i, data_j, i, j);
        p.template add<factor>(res.distance, res.weight);
    }
}