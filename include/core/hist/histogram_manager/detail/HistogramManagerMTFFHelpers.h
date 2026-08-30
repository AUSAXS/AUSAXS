// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelperBase.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>

namespace ausaxs::hist::detail {
    template<bool vbw, int factor>
    inline void evaluate_aa16(WeightedDistribution3D& p, const CompactCoordinatesFF<vbw>& data_a, int i, int j) {
        xyzff::HexaEvaluatedResult res = add16::evaluate_weighted(data_a, data_a, i, j);
        for (int k = 0; k < 16; ++k) {
            p.increment_linear_index<factor>(res.ff_bins[k], res.distance_bins[k], res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aa16(Distribution3D& p, const CompactCoordinatesFF<vbw>& data_a, int i, int j) {
        xyzff::HexaEvaluatedResultRounded res = add16::evaluate_unweighted(data_a, data_a, i, j);
        for (int k = 0; k < 16; ++k) {
            p.increment_linear_index<factor>(res.ff_bins[k], res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aa8(WeightedDistribution3D& p, const CompactCoordinatesFF<vbw>& data_a, int i, int j) {
        xyzff::OctoEvaluatedResult res = add8::evaluate_weighted(data_a, data_a, i, j);
        for (int k = 0; k < 8; ++k) {
            p.increment_linear_index<factor>(res.ff_bins[k], res.distance_bins[k], res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aa8(Distribution3D& p, const CompactCoordinatesFF<vbw>& data_a, int i, int j) {
        xyzff::OctoEvaluatedResultRounded res = add8::evaluate_unweighted(data_a, data_a, i, j);
        for (int k = 0; k < 8; ++k) {
            p.increment_linear_index<factor>(res.ff_bins[k], res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aa4(WeightedDistribution3D& p, const CompactCoordinatesFF<vbw>& data_a, int i, int j) {
        xyzff::QuadEvaluatedResult res = add4::evaluate_weighted(data_a, data_a, i, j);
        for (int k = 0; k < 4; ++k) {
            p.increment_linear_index<factor>(res.ff_bins[k], res.distance_bins[k], res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aa4(Distribution3D& p, const CompactCoordinatesFF<vbw>& data_a, int i, int j) {
        xyzff::QuadEvaluatedResultRounded res = add4::evaluate_unweighted(data_a, data_a, i, j);
        for (int k = 0; k < 4; ++k) {
            p.increment_linear_index<factor>(res.ff_bins[k], res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aa1(WeightedDistribution3D& p, const CompactCoordinatesFF<vbw>& data_a, int i, int j) {
        xyzff::EvaluatedResult res = add1::evaluate_weighted(data_a, data_a, i, j);
        p.increment_linear_index<factor>(res.ff_bin, res.distance_bin, res.distance);
    }

    template<bool vbw, int factor>
    inline void evaluate_aa1(Distribution3D& p, const CompactCoordinatesFF<vbw>& data_a, int i, int j) {
        xyzff::EvaluatedResultRounded res = add1::evaluate_unweighted(data_a, data_a, i, j);
        p.increment_linear_index<factor>(res.ff_bin, res.distance);
    }

    template<bool vbw, int factor>
    inline void evaluate_aw16(WeightedDistribution2D& p, const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j) {
        xyzff::HexaEvaluatedResult res = add16::evaluate_weighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 16; ++k) {
            p.increment_index<factor>(ff_i, res.distance_bins[k], res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aw16(Distribution2D& p, const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j) {
        xyzff::HexaEvaluatedResultRounded res = add16::evaluate_unweighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 16; ++k) {
            p.increment_index<factor>(ff_i, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aw8(WeightedDistribution2D& p, const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j) {
        xyzff::OctoEvaluatedResult res = add8::evaluate_weighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 8; ++k) {
            p.increment_index<factor>(ff_i, res.distance_bins[k], res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aw8(Distribution2D& p, const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j) {
        xyzff::OctoEvaluatedResultRounded res = add8::evaluate_unweighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 8; ++k) {
            p.increment_index<factor>(ff_i, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aw4(WeightedDistribution2D& p, const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j) {
        xyzff::QuadEvaluatedResult res = add4::evaluate_weighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 4; ++k) {
            p.increment_index<factor>(ff_i, res.distance_bins[k], res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aw4(Distribution2D& p, const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j) {
        xyzff::QuadEvaluatedResultRounded res = add4::evaluate_unweighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 4; ++k) {
            p.increment_index<factor>(ff_i, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    inline void evaluate_aw1(WeightedDistribution2D& p, const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j) {
        xyzff::EvaluatedResult res = add1::evaluate_weighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        p.increment_index<factor>(ff_i, res.distance_bin, res.distance);
    }

    template<bool vbw, int factor>
    inline void evaluate_aw1(Distribution2D& p, const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j) {
        xyzff::EvaluatedResultRounded res = add1::evaluate_unweighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        p.increment_index<factor>(ff_i, res.distance);
    }
}