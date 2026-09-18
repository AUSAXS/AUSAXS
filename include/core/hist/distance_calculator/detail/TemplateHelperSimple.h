// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/distance_calculator/detail/TemplateHelperBase.h>
#include <hist/distribution/detail/WeightedEntry.h>

#include <cassert>
#include <cstdint>
#include <span>
#include <type_traits>

namespace ausaxs::hist::detail {
    /**
     * @brief A bin of either kind of distance distribution: a bare value, or one that also tracks the weighted
     *        centre of the bin.
     */
    template<typename Entry>
    concept BinEntry = std::is_same_v<Entry, constants::axes::d_type> || std::is_same_v<Entry, WeightedEntry>;

    // The bins are taken as a span rather than as the distribution that owns them, so that these also serve a
    // single row of a 2D or 3D distribution - which is what a calculation with fixed form factor indices
    // accumulates into. The element type is what selects the interface.

    /**
     * @brief The bins of a whole distribution, as the evaluators take them.
     */
    template<typename Distribution>
    inline std::span<typename Distribution::value_type> bins(Distribution& distribution) {
        return {&*distribution.begin(), static_cast<std::size_t>(distribution.size())};
    }

    template<int factor, BinEntry Entry>
    inline void accumulate(std::span<Entry> p, std::int32_t bin, float distance, float weight) {
        assert(0 <= bin && bin < static_cast<std::int32_t>(p.size()) && "hist::detail::accumulate: bin index out of bounds.");
        if constexpr (std::is_same_v<Entry, WeightedEntry>) {
            p[bin].template add<factor>(distance, weight);
        } else {
            p[bin] += factor*weight;
        }
    }

    template<bool variable_bin_widths, int factor, BinEntry Entry>
    inline void evaluate16(std::span<Entry> p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        if constexpr (std::is_same_v<Entry, WeightedEntry>) {
            xyzw::HexaEvaluatedResult res = add16::evaluate_weighted(data_i, data_j, i, j);
            for (int k = 0; k < 16; ++k) {accumulate<factor>(p, res.distance_bins[k], res.distances[k], res.weights[k]);}
        } else {
            xyzw::HexaEvaluatedResultRounded res = add16::evaluate_unweighted(data_i, data_j, i, j);
            for (int k = 0; k < 16; ++k) {accumulate<factor>(p, res.distance_bins[k], 0, res.weights[k]);}
        }
    }

    template<bool variable_bin_widths, int factor, BinEntry Entry>
    inline void evaluate8(std::span<Entry> p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        if constexpr (std::is_same_v<Entry, WeightedEntry>) {
            xyzw::OctoEvaluatedResult res = add8::evaluate_weighted(data_i, data_j, i, j);
            for (int k = 0; k < 8; ++k) {accumulate<factor>(p, res.distance_bins[k], res.distances[k], res.weights[k]);}
        } else {
            xyzw::OctoEvaluatedResultRounded res = add8::evaluate_unweighted(data_i, data_j, i, j);
            for (int k = 0; k < 8; ++k) {accumulate<factor>(p, res.distance_bins[k], 0, res.weights[k]);}
        }
    }

    template<bool variable_bin_widths, int factor, BinEntry Entry>
    inline void evaluate4(std::span<Entry> p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        if constexpr (std::is_same_v<Entry, WeightedEntry>) {
            xyzw::QuadEvaluatedResult res = add4::evaluate_weighted(data_i, data_j, i, j);
            for (int k = 0; k < 4; ++k) {accumulate<factor>(p, res.distance_bins[k], res.distances[k], res.weights[k]);}
        } else {
            xyzw::QuadEvaluatedResultRounded res = add4::evaluate_unweighted(data_i, data_j, i, j);
            for (int k = 0; k < 4; ++k) {accumulate<factor>(p, res.distance_bins[k], 0, res.weights[k]);}
        }
    }

    template<bool variable_bin_widths, int factor, BinEntry Entry>
    inline void evaluate1(std::span<Entry> p, const CompactCoordinates<variable_bin_widths>& data_i, const CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        if constexpr (std::is_same_v<Entry, WeightedEntry>) {
            xyzw::EvaluatedResult res = add1::evaluate_weighted(data_i, data_j, i, j);
            accumulate<factor>(p, res.distance_bin, res.distance, res.weight);
        } else {
            xyzw::EvaluatedResultRounded res = add1::evaluate_unweighted(data_i, data_j, i, j);
            accumulate<factor>(p, res.distance_bin, 0, res.weight);
        }
    }
}
