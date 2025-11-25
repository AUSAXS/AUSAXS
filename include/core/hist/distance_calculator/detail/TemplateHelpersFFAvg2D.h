// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <form_factor/FormFactorType.h>

namespace ausaxs {
    /**
     * @brief Calculate the distances between eight atoms and template add them to the histogram. The contribution of excluded volume dummy atoms is template added to the excluded volume bin.
     * 
     * @tparam weighted_bins Whether to keep track of the distances template added to the bins. This is useful for weighting the bins later.
     * @tparam factor A multiplicative factor for the atomic weights. 
     * @param p The histogram to template add the distances to.
     * @param data_i The first atom. The form factor index of this will be used for the first axis of the histogram.
     * @param data_j The second atom. This is assumed to be water.
     * @param i The index of the first atom.
     * @param j The index of the second atom.
     */
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution2D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add8::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 8; ++k) {
            if constexpr (factor == 1) {
                p.template increment<1>(res.ff_bins[k], res.distances[k]);
                p.template increment<1>(exv_bin, res.distances[k]);
            } else {
                p.template increment<2>(res.ff_bins[k], res.distances[k]);
                p.template increment<2>(exv_bin, res.distances[k]);
            }
        }
    }

    /**
     * @brief Calculate the distances between four atoms and template add them to the histogram. The contribution of excluded volume dummy atoms is template added to the excluded volume bin.
     * 
     * @tparam weighted_bins Whether to keep track of the distances template added to the bins. This is useful for weighting the bins later.
     * @tparam factor A multiplicative factor for the atomic weights. 
     * @param p The histogram to template add the distances to.
     * @param data_i The first atom. The form factor index of this will be used for the first axis of the histogram.
     * @param data_j The second atom. This is assumed to be water.
     * @param i The index of the first atom.
     * @param j The index of the second atom.
     */
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution2D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add4::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 4; ++k) {
            if constexpr (factor == 1) {
                p.template increment<1>(res.ff_bins[k], res.distances[k]);
                p.template increment<1>(exv_bin, res.distances[k]);
            } else {
                p.template increment<2>(res.ff_bins[k], res.distances[k]);
                p.template increment<2>(exv_bin, res.distances[k]);
            }
        }
    }

    /**
     * @brief Calculate the distances between two atoms and template add them to the histogram. The contribution of excluded volume dummy atoms is template added to the excluded volume bin.
     * 
     * @tparam weighted_bins Whether to keep track of the distances template added to the bins. This is useful for weighting the bins later.
     * @tparam factor A multiplicative factor for the atomic weights. 
     * @param p The histogram to template add the distances to.
     * @param data_i The first atom. The form factor index of this will be used for the first axis of the histogram.
     * @param data_j The second atom. This is assumed to be water.
     * @param i The index of the first atom.
     * @param j The index of the second atom.
     */
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution2D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add1::evaluate<weighted_bins>(data_i, data_j, i, j);
        if constexpr (factor == 1) {
            p.template increment<1>(res.ff_bin, res.distance);
            p.template increment<1>(exv_bin, res.distance);
        } else {
            p.template increment<2>(res.ff_bin, res.distance);
            p.template increment<2>(exv_bin, res.distance);
        }
    }
}