// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>

// Base template helpers - core distance evaluation functions used by all histogram managers
// These provide the low-level SIMD distance calculations
// Separate functions for weighted (float distances) vs unweighted (int32_t bins)

namespace ausaxs::hist::detail {
    namespace add8 {
        //### CompactCoordinates overloads ###//
        template<bool variable_bin_width>
        inline hist::detail::xyzw::OctoEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
        }

        template<bool variable_bin_width>
        inline hist::detail::xyzw::OctoEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool variable_bin_width>
        inline hist::detail::xyzff::OctoEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
        }

        template<bool variable_bin_width>
        inline hist::detail::xyzff::OctoEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
        }
    }

    namespace add4 {
        //### CompactCoordinates overloads ###//
        template<bool variable_bin_width>
        inline hist::detail::xyzw::QuadEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
        }

        template<bool variable_bin_width>
        inline hist::detail::xyzw::QuadEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool variable_bin_width>
        inline hist::detail::xyzff::QuadEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
        }

        template<bool variable_bin_width>
        inline hist::detail::xyzff::QuadEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
        }
    }

    namespace add1 {
        //### CompactCoordinates overloads ###//
        template<bool variable_bin_width>
        inline hist::detail::xyzw::EvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate(data_j[j]);
        }

        template<bool variable_bin_width>
        inline hist::detail::xyzw::EvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded(data_j[j]);
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool variable_bin_width>
        inline hist::detail::xyzff::EvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate(data_j[j]);
        }

        template<bool variable_bin_width>
        inline hist::detail::xyzff::EvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded(data_j[j]);
        }
    }
}

//! The remaining should probably be deleted

// Mixed-type overloads for CompactCoordinatesFF + CompactCoordinates (atom/water + excluded volume)
// For FF-based calculations, we only need distances - weights are irrelevant since FF scaling is applied separately.
// Returns xyzw results with weight=1.0 for each distance.
namespace ausaxs::hist::detail::add8 {
    template<bool variable_bin_width>
    inline auto evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        // Create temporary CompactCoordinates from the XYZFF's position with weight 1.0
        hist::detail::CompactCoordinatesXYZW<variable_bin_width> temp_i(data_i[i].value.pos, 1.0f);
        return temp_i.evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
    }

    template<bool variable_bin_width>
    inline auto evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        hist::detail::CompactCoordinatesXYZW<variable_bin_width> temp_i(data_i[i].value.pos, 1.0f);
        return temp_i.evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3], data_j[j+4], data_j[j+5], data_j[j+6], data_j[j+7]);
    }
}

namespace ausaxs::hist::detail::add4 {
    template<bool variable_bin_width>
    inline auto evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        hist::detail::CompactCoordinatesXYZW<variable_bin_width> temp_i(data_i[i].value.pos, 1.0f);
        return temp_i.evaluate(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
    }

    template<bool variable_bin_width>
    inline auto evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        hist::detail::CompactCoordinatesXYZW<variable_bin_width> temp_i(data_i[i].value.pos, 1.0f);
        return temp_i.evaluate_rounded(data_j[j], data_j[j+1], data_j[j+2], data_j[j+3]);
    }
}

namespace ausaxs::hist::detail::add1 {
    template<bool variable_bin_width>
    inline auto evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        hist::detail::CompactCoordinatesXYZW<variable_bin_width> temp_i(data_i[i].value.pos, 1.0f);
        return temp_i.evaluate(data_j[j]);
    }

    template<bool variable_bin_width>
    inline auto evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
        hist::detail::CompactCoordinatesXYZW<variable_bin_width> temp_i(data_i[i].value.pos, 1.0f);
        return temp_i.evaluate_rounded(data_j[j]);
    }
}
