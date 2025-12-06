// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelperBase.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distribution/GenericDistribution1D.h>

// Template helpers for basic 1D histogram managers
// Used by: HistogramManager, HistogramManagerMT, PartialHistogramManager, 
//          PartialHistogramManagerMT, PartialSymmetryManagerMT, SimpleCalculator

namespace ausaxs {
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinates<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add8::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 8; ++k) {
            p.template add<factor>(res.distances[k], res.weights[k]);
        }
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinates<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add4::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 4; ++k) {
            p.template add<factor>(res.distances[k], res.weights[k]);
        }
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinates<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add1::evaluate<weighted_bins>(data_i, data_j, i, j);
        p.template add<factor>(res.distance, res.weight);
    }
}
