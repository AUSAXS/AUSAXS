// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distance_calculator/detail/TemplateHelpersFFAvg2D.h>
#include <hist/distance_calculator/detail/TemplateHelpersFFAvg3D.h>

// Mixed-type helpers for CompactCoordinatesFF with CompactCoordinates (water with excluded volume)
// These simply delegate to the base CompactCoordinates evaluate methods since water's form factor is irrelevant here
namespace ausaxs {
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        // Reinterpret CompactCoordinatesFF as CompactCoordinates (same memory layout, different template parameter)
        evaluate8<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), data_j, i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        evaluate4<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), data_j, i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinates<variable_bin_widths>& data_j, int i, int j) {
        evaluate1<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), data_j, i, j);
    }

    // Water-water helpers (CompactCoordinatesFF + CompactCoordinatesFF → 1D histogram)
    // Water has constant form factor, so we just compute distances
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        // Reinterpret both as CompactCoordinates since we only care about distances
        evaluate8<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_j), i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        evaluate4<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_j), i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        evaluate1<weighted_bins, variable_bin_widths, factor>(p, reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_j), i, j);
    }
}