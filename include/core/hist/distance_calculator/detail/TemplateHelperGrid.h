// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelperBase.h>
#include <hist/distance_calculator/detail/TemplateHelperAvg.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/detail/CompactCoordinatesFF.h>

// Template helpers for FFGrid histogram managers
// Used by: HistogramManagerMTFFGrid, HistogramManagerMTFFGridSurface, HistogramManagerMTFFGridScalableExv

// Water-water overloads for 1D histograms (CompactCoordinatesFF + CompactCoordinatesFF → 1D)
// Water molecules all have the same form factor type, so we just need distances.
// For weighted bins, all waters have effective weight 1.0.
namespace ausaxs {
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add8::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 8; ++k) {
            if constexpr (weighted_bins) {p.template add<factor>(res.distances[k], 1.0f);}
            else {p.template increment<factor>(res.distances[k]);}
        }
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add4::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 4; ++k) {
            if constexpr (weighted_bins) {p.template add<factor>(res.distances[k], 1.0f);}
            else {p.template increment<factor>(res.distances[k]);}
        }
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add1::evaluate<weighted_bins>(data_i, data_j, i, j);
        if constexpr (weighted_bins) {p.template add<factor>(res.distance, 1.0f);}
        else {p.template increment<factor>(res.distance);}
    }

    // Water-excluded volume overloads (CompactCoordinatesFF + CompactCoordinates → 1D)
    // Water has constant form factor, so just compute weighted distances
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add8::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 8; ++k) {
            if constexpr (weighted_bins) {p.template add<factor>(res.distances[k], res.weights[k]);}
            else {p.template increment<factor>(res.distances[k]);}
        }
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add4::evaluate<weighted_bins>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 4; ++k) {
            if constexpr (weighted_bins) {p.template add<factor>(res.distances[k], res.weights[k]);}
            else {p.template increment<factor>(res.distances[k]);}
        }
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = detail::add1::evaluate<weighted_bins>(data_i, data_j, i, j);
        if constexpr (weighted_bins) {p.template add<factor>(res.distance, res.weight);}
        else {p.template increment<factor>(res.distance);}
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
