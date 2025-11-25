// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelperBase.h>
#include <hist/distance_calculator/detail/TemplateHelperAvg.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/detail/CompactCoordinatesFF.h>

// Template helpers for FFGrid histogram managers
// Used by: HistogramManagerMTFFGrid, HistogramManagerMTFFGridSurface, HistogramManagerMTFFGridScalableExv

// Mixed-type overloads for water-excluded volume (1D histograms)
// Water has constant form factor, so we just need distances (reinterpret FF as base CompactCoordinates)
namespace ausaxs {
    // CompactCoordinatesFF + CompactCoordinatesFF → 1D
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        evaluate8<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_j), i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        evaluate4<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_j), i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        evaluate1<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_j), i, j);
    }

    // CompactCoordinatesFF + CompactCoordinates → 1D
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        evaluate8<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), data_j, i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        evaluate4<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), data_j, i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        evaluate1<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), data_j, i, j);
    }
}

// Grid-specific evaluates for atom-excluded volume (2D: ff, distance)
// These DO NOT add excluded volume contributions (unlike FFAvg helpers)
namespace ausaxs::grid {
    template<bool weighted_bins, bool variable_bin_width, int factor>
    inline void evaluate8(
        typename hist::GenericDistribution2D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, 
        const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j
    ) {
        auto res = ausaxs::detail::add8::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 8; ++k) {p.add(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool weighted_bins, bool variable_bin_width, int factor>
    inline void evaluate4(
        typename hist::GenericDistribution2D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, 
        const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j
    ) {
        auto res = ausaxs::detail::add4::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 4; ++k) {p.add(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool weighted_bins, bool variable_bin_width, int factor>
    inline void evaluate1(
        typename hist::GenericDistribution2D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, 
        const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j
    ) {
        auto res = ausaxs::detail::add1::evaluate<weighted_bins>(data_i, data_j, i, j);
        p.add(data_i.get_ff_type(i), res.distance, factor*res.weight);
    }
}
