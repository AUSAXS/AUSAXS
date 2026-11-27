// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

// the lattice transform is only available when built with pocketfft; see the POCKETFFT option
#if defined(POCKETFFT_AVAILABLE)

#include <grid/detail/GridExcludedVolume.h>
#include <hist/distribution/WeightedDistribution1D.h>

#include <array>

/**
 * @brief Pair-distance histograms of point sets supported on an orthogonal lattice.
 */
namespace ausaxs::hist::detail::lattice {
    /**
     * @brief The self- and cross-correlations of a two-component lattice point set.
     */
    struct Correlations {
        WeightedDistribution1D first;   // first-first pairs
        WeightedDistribution1D second;  // second-second pairs
        WeightedDistribution1D cross;   // first-second pairs, counted once in each direction
    };

    /**
     * @brief The self-correlation histogram of the interior excluded volume points.
     *
     * @param exv The excluded volume. Only its interior sites and the lattice spacing are used.
     * @param inv_bin_width The inverse histogram bin width, as used by the distance calculators.
     * @param bin_count The size of the returned histogram.
     */
    WeightedDistribution1D self_correlation(const grid::exv::GridExcludedVolume& exv, double inv_bin_width, int bin_count);

    /**
     * @brief The self-correlation histogram of lattice sites with independent axis spacings.
     *
     * Sites must be unique and non-negative. Zero-displacement self-pairs are excluded.
     * @param sites Integer lattice coordinates.
     * @param spacing Positive lattice spacings along x, y and z, in Å.
     * @param inv_bin_width The inverse histogram bin width.
     * @param bin_count The size of the returned histogram.
     */
    WeightedDistribution1D self_correlation(
        const std::vector<Vector3<int>>& sites, const std::array<double, 3>& spacing, double inv_bin_width, int bin_count
    );

    /**
     * @brief The three correlation histograms of the interior and surface excluded volume points.
     *
     * @param exv The excluded volume. Its interior and surface sites and the lattice spacing are used.
     * @param inv_bin_width The inverse histogram bin width, as used by the distance calculators.
     * @param bin_count The size of the returned histograms.
     */
    Correlations correlations(const grid::exv::GridExcludedVolume& exv, double inv_bin_width, int bin_count);
}

#endif
