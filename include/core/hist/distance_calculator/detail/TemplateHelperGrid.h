// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelperBase.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/detail/CompactCoordinatesFF.h>

namespace ausaxs::hist::detail {
    // Weighted distribution overloads
    template<bool variable_bin_widths, bool explicit_ff = false, int factor = 1>
    inline void evaluate8(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths, explicit_ff>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzff::OctoEvaluatedResult res = add8::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.increment_bin<factor>(res.distance_bins[k], res.distances[k]);}
    }

    template<bool variable_bin_widths, bool explicit_ff = false, int factor = 1>
    inline void evaluate4(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths, explicit_ff>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzff::QuadEvaluatedResult res = add4::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.increment_bin<factor>(res.distance_bins[k], res.distances[k]);}
    }

    template<bool variable_bin_widths, bool explicit_ff = false, int factor = 1>
    inline void evaluate1(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths, explicit_ff>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzff::EvaluatedResult res = add1::evaluate_weighted(data_i, data_j, i, j);
        p.increment_bin<factor>(res.distance_bin, res.distance);
    }

    // Unweighted distribution overloads
    template<bool variable_bin_widths, bool explicit_ff = false, int factor = 1>
    inline void evaluate8(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths, explicit_ff>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzff::OctoEvaluatedResultRounded res = add8::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.increment_bin<factor>(res.distances[k]);}
    }

    template<bool variable_bin_widths, bool explicit_ff = false, int factor = 1>
    inline void evaluate4(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths, explicit_ff>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzff::QuadEvaluatedResultRounded res = add4::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.template increment_bin<factor>(res.distances[k]);}
    }

    template<bool variable_bin_widths, bool explicit_ff = false, int factor = 1>
    inline void evaluate1(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths, explicit_ff>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzff::EvaluatedResultRounded res = add1::evaluate_unweighted(data_i, data_j, i, j);
        p.template increment_bin<factor>(res.distance);
    }
}

namespace ausaxs::hist::detail::grid {
    // Weighted distribution overloads - receive float distances
    template<bool variable_bin_width, bool explicit_ff = false, int factor = 1>
    inline void evaluate8(
        WeightedDistribution2D& p, const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, 
        const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j
    ) {
        xyzff::OctoEvaluatedResult res = add8::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.increment_bin<factor>(data_i.get_ff_type(i), res.distance_bins[k], res.distances[k]);}
    }

    template<bool variable_bin_width, bool explicit_ff = false, int factor = 1>
    inline void evaluate4(
        WeightedDistribution2D& p, const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, 
        const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j
    ) {
        xyzff::QuadEvaluatedResult res = add4::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.increment_bin<factor>(data_i.get_ff_type(i), res.distance_bins[k], res.distances[k]);}
    }

    template<bool variable_bin_width, bool explicit_ff = false, int factor = 1>
    inline void evaluate1(
        WeightedDistribution2D& p, const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, 
        const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j
    ) {
        xyzff::EvaluatedResult res = add1::evaluate_weighted(data_i, data_j, i, j);
        p.increment_bin<factor>(data_i.get_ff_type(i), res.distance_bin, res.distance);
    }

    // Unweighted distribution overloads - receive int32_t bins
    template<bool variable_bin_width, bool explicit_ff = false, int factor = 1>
    inline void evaluate8(
        Distribution2D& p, const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, 
        const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j
    ) {
        xyzff::OctoEvaluatedResultRounded res = add8::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.increment_bin<factor>(data_i.get_ff_type(i), res.distances[k]);}
    }

    template<bool variable_bin_width, bool explicit_ff = false, int factor = 1>
    inline void evaluate4(
        Distribution2D& p, const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, 
        const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j
    ) {
        xyzff::QuadEvaluatedResultRounded res = add4::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.increment_bin<factor>(data_i.get_ff_type(i), res.distances[k]);}
    }

    template<bool variable_bin_width, bool explicit_ff = false, int factor = 1>
    inline void evaluate1(
        Distribution2D& p, const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, 
        const CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j
    ) {
        xyzff::EvaluatedResultRounded res = add1::evaluate_unweighted(data_i, data_j, i, j);
        p.increment_bin<factor>(data_i.get_ff_type(i), res.distance);
    }
}
