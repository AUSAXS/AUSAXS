// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelperBase.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace ausaxs::hist::detail {
    //### Weighted evaluators ###//
    template<bool variable_bin_widths, int factor>
    inline void evaluate8(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        xyzff::OctoEvaluatedResult res = add8::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.increment<factor>(res.distances[k]);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate4(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        xyzff::QuadEvaluatedResult res = add4::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.increment<factor>(res.distances[k]);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate1(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        xyzff::EvaluatedResult res = add1::evaluate_weighted(data_i, data_j, i, j);
        p.increment<factor>(res.distance);
    }

    //### Unweighted evaluators ###//
    template<bool variable_bin_widths, int factor>
    inline void evaluate8(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        xyzff::OctoEvaluatedResultRounded res = add8::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.increment_bin<factor>(res.distances[k]);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate4(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        xyzff::QuadEvaluatedResultRounded res = add4::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {
            p.increment_bin<factor>(res.distances[k]);
        }
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate1(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        xyzff::EvaluatedResultRounded res = add1::evaluate_unweighted(data_i, data_j, i, j);
        p.increment_bin<factor>(res.distance);
    }
}
