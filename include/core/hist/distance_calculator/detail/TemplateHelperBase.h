// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>

// This is the seam between the pair loops, which think in (container, index) pairs, and the
// kernels, which take one broadcast atom plus one pointer per component of the opposing block.
// Nothing above this file needs to know how the coordinates are laid out.
namespace ausaxs::hist::detail {
    namespace add16 {
        //### CompactCoordinates overloads ###//
        template<bool vbw>
        inline xyzw::HexaEvaluatedResult evaluate_weighted(const CompactCoordinates<vbw>& data_i, const CompactCoordinates<vbw>& data_j, int i, int j) {
            return xyzw::evaluate_16<vbw>(data_i.atom(i), data_j.block(j));
        }

        template<bool vbw>
        inline xyzw::HexaEvaluatedResultRounded evaluate_unweighted(const CompactCoordinates<vbw>& data_i, const CompactCoordinates<vbw>& data_j, int i, int j) {
            return xyzw::evaluate_rounded_16<vbw>(data_i.atom(i), data_j.block(j));
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool vbw>
        inline xyzff::HexaEvaluatedResult evaluate_weighted(const CompactCoordinatesFF<vbw>& data_i, const CompactCoordinatesFF<vbw>& data_j, int i, int j) {
            return xyzff::evaluate_16<vbw>(data_i.atom(i), data_j.block(j));
        }

        template<bool vbw>
        inline xyzff::HexaEvaluatedResultRounded evaluate_unweighted(const CompactCoordinatesFF<vbw>& data_i, const CompactCoordinatesFF<vbw>& data_j, int i, int j) {
            return xyzff::evaluate_rounded_16<vbw>(data_i.atom(i), data_j.block(j));
        }
    }

    namespace add8 {
        //### CompactCoordinates overloads ###//
        template<bool vbw>
        inline xyzw::OctoEvaluatedResult evaluate_weighted(const CompactCoordinates<vbw>& data_i, const CompactCoordinates<vbw>& data_j, int i, int j) {
            return xyzw::evaluate_8<vbw>(data_i.atom(i), data_j.block(j));
        }

        template<bool vbw>
        inline xyzw::OctoEvaluatedResultRounded evaluate_unweighted(const CompactCoordinates<vbw>& data_i, const CompactCoordinates<vbw>& data_j, int i, int j) {
            return xyzw::evaluate_rounded_8<vbw>(data_i.atom(i), data_j.block(j));
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool vbw>
        inline xyzff::OctoEvaluatedResult evaluate_weighted(const CompactCoordinatesFF<vbw>& data_i, const CompactCoordinatesFF<vbw>& data_j, int i, int j) {
            return xyzff::evaluate_8<vbw>(data_i.atom(i), data_j.block(j));
        }

        template<bool vbw>
        inline xyzff::OctoEvaluatedResultRounded evaluate_unweighted(const CompactCoordinatesFF<vbw>& data_i, const CompactCoordinatesFF<vbw>& data_j, int i, int j) {
            return xyzff::evaluate_rounded_8<vbw>(data_i.atom(i), data_j.block(j));
        }
    }

    namespace add4 {
        //### CompactCoordinates overloads ###//
        template<bool vbw>
        inline xyzw::QuadEvaluatedResult evaluate_weighted(const CompactCoordinates<vbw>& data_i, const CompactCoordinates<vbw>& data_j, int i, int j) {
            return xyzw::evaluate_4<vbw>(data_i.atom(i), data_j.block(j));
        }

        template<bool vbw>
        inline xyzw::QuadEvaluatedResultRounded evaluate_unweighted(const CompactCoordinates<vbw>& data_i, const CompactCoordinates<vbw>& data_j, int i, int j) {
            return xyzw::evaluate_rounded_4<vbw>(data_i.atom(i), data_j.block(j));
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool vbw>
        inline xyzff::QuadEvaluatedResult evaluate_weighted(const CompactCoordinatesFF<vbw>& data_i, const CompactCoordinatesFF<vbw>& data_j, int i, int j) {
            return xyzff::evaluate_4<vbw>(data_i.atom(i), data_j.block(j));
        }

        template<bool vbw>
        inline xyzff::QuadEvaluatedResultRounded evaluate_unweighted(const CompactCoordinatesFF<vbw>& data_i, const CompactCoordinatesFF<vbw>& data_j, int i, int j) {
            return xyzff::evaluate_rounded_4<vbw>(data_i.atom(i), data_j.block(j));
        }
    }

    namespace add1 {
        //### CompactCoordinates overloads ###//
        template<bool vbw>
        inline xyzw::EvaluatedResult evaluate_weighted(const CompactCoordinates<vbw>& data_i, const CompactCoordinates<vbw>& data_j, int i, int j) {
            return xyzw::evaluate<vbw>(data_i.atom(i), data_j.block(j));
        }

        template<bool vbw>
        inline xyzw::EvaluatedResultRounded evaluate_unweighted(const CompactCoordinates<vbw>& data_i, const CompactCoordinates<vbw>& data_j, int i, int j) {
            return xyzw::evaluate_rounded<vbw>(data_i.atom(i), data_j.block(j));
        }

        //### CompactCoordinatesFF overloads ###//
        template<bool vbw>
        inline xyzff::EvaluatedResult evaluate_weighted(const CompactCoordinatesFF<vbw>& data_i, const CompactCoordinatesFF<vbw>& data_j, int i, int j) {
            return xyzff::evaluate<vbw>(data_i.atom(i), data_j.block(j));
        }

        template<bool vbw>
        inline xyzff::EvaluatedResultRounded evaluate_unweighted(const CompactCoordinatesFF<vbw>& data_i, const CompactCoordinatesFF<vbw>& data_j, int i, int j) {
            return xyzff::evaluate_rounded<vbw>(data_i.atom(i), data_j.block(j));
        }
    }
}
