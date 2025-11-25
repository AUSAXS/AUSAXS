// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distance_calculator/detail/TemplateHelpersFFExplicit2D.h>
#include <hist/distance_calculator/detail/TemplateHelpersFFExplicit3D.h>

// Water-water helpers (CompactCoordinatesFF + CompactCoordinatesFF → 1D histogram)
// Water has constant form factor, so we just compute distances
namespace ausaxs {
    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate8(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        // Cast both to base since we only care about distances
        evaluate8<weighted_bins, variable_bin_widths, factor>(p, static_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), static_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_j), i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate4(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        evaluate4<weighted_bins, variable_bin_widths, factor>(p, static_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), static_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_j), i, j);
    }

    template<bool weighted_bins, bool variable_bin_widths, int factor>
    inline void evaluate1(typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_i, const hist::detail::CompactCoordinatesFF<variable_bin_widths>& data_j, int i, int j) {
        evaluate1<weighted_bins, variable_bin_widths, factor>(p, static_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_i), static_cast<const hist::detail::CompactCoordinates<variable_bin_widths>&>(data_j), i, j);
    }
}