// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/WeightedDistribution1D.h>
#include <math/MathFwd.h>

#include <optional>
#include <vector>

/**
 * @brief Pair-distance histograms of point sets supported on a cubic lattice.
 *
 * The grid-based excluded volume models emit voxel centers, so their points are not arbitrary: they are an exact
 * subset of the sites of a single cubic lattice. The self-correlation of such a set is the radially binned
 * autocorrelation of a binary occupancy array, which the Wiener-Khinchin theorem turns into three transforms of the
 * zero-padded bounding box. That is O(M log M) in the padded box volume M in place of O(N^2) in the point count, and
 * since the excluded volume holds ~17 points per atom, the pair count it replaces is ~289x that of the atom-atom loop.
 *
 * The result is *exact*, not an approximation: zero-padding each axis to at least 2*extent-1 removes the circular
 * wrap-around, every displacement has an exactly known integer squared length, and the counts come out as integers.
 * The distances are in fact more accurate than those of the pair loop, which works in float32 throughout.
 *
 * The cost is memory - the transform buffer scales with the *bounding box* rather than the point count, so it grows
 * with the bounding volume rather than the occupied one. That buffer is a single copy of the padded box, 8 bytes per
 * cell: both transforms run in place, and the several spellings that silently cost two or three copies of it are
 * avoided deliberately. See the Transform class comment in the implementation for how, and why the obvious float32
 * halving is a trap rather than a further 2x. At one buffer the requirement is small enough for ordinary structures
 * that it is simply allocated, exactly as grid::Grid allocates what its own axes need.
 */
namespace ausaxs::hist::detail::lattice {
    /**
     * @brief The self- and cross-correlations of a two-component lattice point set.
     *
     * All three are ordered pair counts, matching what the equivalent pair loops using evaluate<..., 2> produce.
     * Self-pairs are excluded from @a first and @a second, as the callers add those terms themselves.
     */
    struct Correlations {
        WeightedDistribution1D first;   // first-first pairs
        WeightedDistribution1D second;  // second-second pairs
        WeightedDistribution1D cross;   // first-second pairs, counted once in each direction
    };

    /**
     * @brief The self-correlation histogram of a point set lying on a cubic lattice of the given spacing.
     *
     * The returned counts are for *ordered* pairs, matching a pair loop using evaluate<..., 2>. Self-pairs are
     * excluded, since the callers add that term themselves.
     *
     * @param points The point set. Must lie on a common cubic lattice of the given spacing.
     * @param spacing The lattice spacing in Ångström. A non-positive value is interpreted as "not lattice-supported".
     * @param inv_bin_width The inverse histogram bin width, as used by the distance calculators.
     * @param bin_count The size of the returned histogram.
     *
     * @return std::nullopt if the points do not lie on the given lattice, or if the transform buffer could not be
     *         allocated. The caller must then fall back to an explicit pair loop.
     */
    std::optional<WeightedDistribution1D> self_correlation(
        const std::vector<Vector3<double>>& points,
        double spacing,
        double inv_bin_width,
        int bin_count
    );

    /**
     * @brief The three correlation histograms of two point sets sharing a cubic lattice of the given spacing.
     *
     * @param first The first point set. Must lie on a common cubic lattice of the given spacing.
     * @param second The second point set. Must lie on the same lattice as @a first.
     * @param spacing The lattice spacing in Ångström. A non-positive value is interpreted as "not lattice-supported".
     * @param inv_bin_width The inverse histogram bin width, as used by the distance calculators.
     * @param bin_count The size of the returned histograms.
     *
     * @return std::nullopt under the same conditions as self_correlation.
     */
    std::optional<Correlations> correlations(
        const std::vector<Vector3<double>>& first,
        const std::vector<Vector3<double>>& second,
        double spacing,
        double inv_bin_width,
        int bin_count
    );
}
