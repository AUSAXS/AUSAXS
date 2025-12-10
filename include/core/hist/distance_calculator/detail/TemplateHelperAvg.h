// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelperBase.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace ausaxs::hist::detail {
    // Weighted distribution overloads: res.distances are exact float distances
    template<bool variable_bin_widths, int factor>
    inline void evaluate8(WeightedDistribution1D& p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzw::OctoEvaluatedResult res = add8::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.add<factor>(res.distances[k], res.weights[k]);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate4(WeightedDistribution1D& p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzw::QuadEvaluatedResult res = add4::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.add<factor>(res.distances[k], res.weights[k]);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate1(WeightedDistribution1D& p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzw::EvaluatedResult res = add1::evaluate_weighted(data_i, data_j, i, j);
        p.add<factor>(res.distance, res.weight);
    }

    // Unweighted distribution overloads: res.distances are int32_t bins
    template<bool variable_bin_widths, int factor>
    inline void evaluate8(Distribution1D& p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzw::OctoEvaluatedResultRounded res = add8::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.add_index<factor>(res.distances[k], res.weights[k]);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate4(Distribution1D& p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzw::QuadEvaluatedResultRounded res = add4::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.add_index<factor>(res.distances[k], res.weights[k]);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate1(Distribution1D& p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        xyzw::EvaluatedResultRounded res = add1::evaluate_unweighted(data_i, data_j, i, j);
        p.add_index<factor>(res.distance, res.weight);
    }
}
