// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>

#include <span>
#include <type_traits>

namespace ausaxs::hist::detail {
    namespace add16 {
        //### CompactCoordinates overloads ###//
        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzw::HexaEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_16(std::span<std::remove_reference_t<decltype(data_j[j])>, 16>(&data_j[j], 16));
        }

        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzw::HexaEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded_16(std::span<std::remove_reference_t<decltype(data_j[j])>, 16>(&data_j[j], 16));
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzff::HexaEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j) {
            return data_i[i].evaluate_16(std::span<std::remove_reference_t<decltype(data_j[j])>, 16>(&data_j[j], 16));
        }

        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzff::HexaEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded_16(std::span<std::remove_reference_t<decltype(data_j[j])>, 16>(&data_j[j], 16));
        }
    }

    namespace add8 {
        //### CompactCoordinates overloads ###//
        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzw::OctoEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_8(std::span<std::remove_reference_t<decltype(data_j[j])>, 8>(&data_j[j], 8));
        }

        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzw::OctoEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded_8(std::span<std::remove_reference_t<decltype(data_j[j])>, 8>(&data_j[j], 8));
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzff::OctoEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j) {
            return data_i[i].evaluate_8(std::span<std::remove_reference_t<decltype(data_j[j])>, 8>(&data_j[j], 8));
        }

        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzff::OctoEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded_8(std::span<std::remove_reference_t<decltype(data_j[j])>, 8>(&data_j[j], 8));
        }
    }

    namespace add4 {
        //### CompactCoordinates overloads ###//
        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzw::QuadEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_4(std::span<std::remove_reference_t<decltype(data_j[j])>, 4>(&data_j[j], 4));
        }

        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzw::QuadEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded_4(std::span<std::remove_reference_t<decltype(data_j[j])>, 4>(&data_j[j], 4));
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzff::QuadEvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j) {
            return data_i[i].evaluate_4(std::span<std::remove_reference_t<decltype(data_j[j])>, 4>(&data_j[j], 4));
        }

        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzff::QuadEvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded_4(std::span<std::remove_reference_t<decltype(data_j[j])>, 4>(&data_j[j], 4));
        }
    }

    namespace add1 {
        //### CompactCoordinates overloads ###//
        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzw::EvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate(data_j[j]);
        }

        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzw::EvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinates<variable_bin_width>& data_i, const hist::detail::CompactCoordinates<variable_bin_width>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded(data_j[j]);
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzff::EvaluatedResult evaluate_weighted(const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j) {
            return data_i[i].evaluate(data_j[j]);
        }

        template<bool variable_bin_width, bool explicit_ff = false>
        inline hist::detail::xyzff::EvaluatedResultRounded evaluate_unweighted(const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_width, explicit_ff>& data_j, int i, int j) {
            return data_i[i].evaluate_rounded(data_j[j]);
        }
    }
}