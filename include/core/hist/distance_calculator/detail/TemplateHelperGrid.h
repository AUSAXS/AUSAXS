// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelperBase.h>
#include <hist/distance_calculator/detail/TemplateHelperAvg.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/detail/CompactCoordinatesFF.h>

namespace ausaxs::hist::detail {
    template<bool variable_bin_widths, int factor>
    inline void evaluate8(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = add8::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.add<factor>(res.distances[k], 1.0f);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate4(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = add4::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.add<factor>(res.distances[k], 1.0f);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate1(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = add1::evaluate_weighted(data_i, data_j, i, j);
        p.add<factor>(res.distance, 1.0f);
    }

    // Unweighted distribution - receives int32_t bins, increments directly
    template<bool variable_bin_widths, int factor>
    inline void evaluate8(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = add8::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.template increment<factor>(res.distances[k]);}
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate4(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = add4::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {
            p.template increment<factor>(res.distances[k]);
        }
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate1(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        auto res = add1::evaluate_unweighted(data_i, data_j, i, j);
        p.template increment<factor>(res.distance);
    }

    // Water-excluded volume overloads (CompactCoordinatesFF + CompactCoordinates → 1D)
    // Water has constant form factor, so compute weighted distances

    // Weighted distribution overloads
    template<bool variable_bin_widths, int factor>
    inline void evaluate8(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = add8::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {
            p.template add<factor>(res.distances[k], res.weights[k]);
        }
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate4(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = add4::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {
            p.template add<factor>(res.distances[k], res.weights[k]);
        }
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate1(WeightedDistribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = add1::evaluate_weighted(data_i, data_j, i, j);
        p.template add<factor>(res.distance, res.weight);
    }

    // Unweighted distribution overloads
    template<bool variable_bin_widths, int factor>
    inline void evaluate8(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = add8::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {
            p.template increment<factor>(res.distances[k]);
        }
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate4(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = add4::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {
            p.template increment<factor>(res.distances[k]);
        }
    }

    template<bool variable_bin_widths, int factor>
    inline void evaluate1(Distribution1D& p, const CompactCoordinatesFF<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        auto res = add1::evaluate_unweighted(data_i, data_j, i, j);
        p.template increment<factor>(res.distance);
    }
}

namespace ausaxs::hist::detail::grid {
    // Weighted distribution overloads - receive float distances
    template<bool variable_bin_width, int factor>
    inline void evaluate8(
        WeightedDistribution2D& p, const CompactCoordinatesFF<variable_bin_width>& data_i, 
        const CompactCoordinates<variable_bin_width>& data_j, int i, int j
    ) {
        auto res = add8::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.add(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool variable_bin_width, int factor>
    inline void evaluate4(
        WeightedDistribution2D& p, const CompactCoordinatesFF<variable_bin_width>& data_i, 
        const CompactCoordinates<variable_bin_width>& data_j, int i, int j
    ) {
        auto res = add4::evaluate_weighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.add(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool variable_bin_width, int factor>
    inline void evaluate1(
        WeightedDistribution2D& p, const CompactCoordinatesFF<variable_bin_width>& data_i, 
        const CompactCoordinates<variable_bin_width>& data_j, int i, int j
    ) {
        auto res = add1::evaluate_weighted(data_i, data_j, i, j);
        p.add(data_i.get_ff_type(i), res.distance, factor*res.weight);
    }

    // Unweighted distribution overloads - receive int32_t bins
    template<bool variable_bin_width, int factor>
    inline void evaluate8(
        Distribution2D& p, const CompactCoordinatesFF<variable_bin_width>& data_i, 
        const CompactCoordinates<variable_bin_width>& data_j, int i, int j
    ) {
        auto res = add8::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 8; ++k) {p.add_index(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool variable_bin_width, int factor>
    inline void evaluate4(
        Distribution2D& p, const CompactCoordinatesFF<variable_bin_width>& data_i, 
        const CompactCoordinates<variable_bin_width>& data_j, int i, int j
    ) {
        auto res = add4::evaluate_unweighted(data_i, data_j, i, j);
        for (int k = 0; k < 4; ++k) {p.add_index(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool variable_bin_width, int factor>
    inline void evaluate1(
        Distribution2D& p, const CompactCoordinatesFF<variable_bin_width>& data_i, 
        const CompactCoordinates<variable_bin_width>& data_j, int i, int j
    ) {
        auto res = add1::evaluate_unweighted(data_i, data_j, i, j);
        p.add_index(data_i.get_ff_type(i), res.distance, factor*res.weight);
    }
}
