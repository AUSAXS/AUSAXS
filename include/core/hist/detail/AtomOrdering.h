// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/detail/WeightedEntry.h>
#include <constants/ConstantsAxes.h>

#include <cstddef>
#include <ranges>

namespace ausaxs::hist::detail {
    namespace atom_order {
        /**
         * @brief The largest histogram footprint for which decorrelating the atom order pays off.
         *
         * Two competing effects decide this, both measured on the accumulation, which
         * dominates the distance-histogram loop:
         *
         * - Linear file order puts spatially adjacent atoms adjacent in memory, so
         *   consecutive inner-loop atoms land in the *same* bin as each other. Each repeat
         *   serialises the read-modify-write on store-to-load forwarding, which is worth
         *   up to 9% of the whole loop. Permuting the order removes it.
         * - That same correlation keeps the set of bins touched at any moment inside a
         *   narrow window, which stays in L1 even when the whole histogram does not.
         *   Permuting the order spreads every block over the entire histogram, so once the
         *   histogram no longer fits in L1 the permutation costs more than it saves.
         *
         * Measured crossover (Zen 5, 4.7k-atom structure, footprint varied by binning the
         * same distances more finely, which is the crystal regime):
         *
         *      footprint   3 kB    11 kB   22 kB   44 kB   88 kB   176 kB
         *      unweighted  -9.1%   -6.2%   -3.0%   -1.5%   -1.5%   +1.1%
         *      weighted    -5.7%   -2.7%   -0.8%   +0.4%   +0.7%   +0.7%
         *
         * 24 kB keeps the whole useful range and leaves a wide margin to the region where
         * it turns into a loss. It is deliberately below a single L1 (32-48 kB on current
         * parts) because two SMT threads share one L1 and the coordinate stream competes
         * for it as well.
         */
        constexpr std::size_t max_histogram_bytes = 24*1024;

        /**
         * @brief The bytes one histogram of @a bin_count bins occupies.
         *        Weighted bins track a running centre and count per bin, so their entries are wider.
         */
        template<bool weighted_bins>
        constexpr std::size_t histogram_bytes(unsigned int bin_count) {
            if constexpr (weighted_bins) {
                return bin_count*sizeof(hist::detail::WeightedEntry);
            } else {
                return bin_count*sizeof(constants::axes::d_type);
            }
        }

        /**
         * @brief Whether decorrelating the atom order is expected to be a net win.
         *
         * Keyed on the histogram footprint rather than on the molecular size, because the
         * bin count is a product of size *and* bin width - a crystal at 0.1 A bins reaches
         * a large histogram from a small structure.
         */
        template<bool weighted_bins>
        constexpr bool is_beneficial(unsigned int bin_count) {
            return histogram_bytes<weighted_bins>(bin_count) <= max_histogram_bytes;
        }

        // a coordinate set that can be permuted directly
        template<typename T>
        concept Shufflable = requires(T& t) {t.shuffle_order(0u);};

        template<Shufflable T>
        void shuffle_all(T& set, unsigned int& seed) {set.shuffle_order(seed++);}

        // a container of sets - possibly nested, as in the per-body symmetry data
        template<std::ranges::input_range Range> requires (!Shufflable<Range>)
        void shuffle_all(Range& sets, unsigned int& seed) {
            for (auto& set : sets) {shuffle_all(set, seed);}
        }
    }

    /**
     * @brief Permute the atom order of the given coordinate sets, if the histogram is small
     *        enough for that to speed up the accumulation.
     *
     * The distance histogram is a sum over all pairs, so this never changes the result
     * beyond floating-point summation order. Each set is permuted independently, which is
     * safe because nothing maps an index in a coordinate set back to the atom it came
     * from - only aggregate quantities are taken from the source bodies.
     *
     * Must be called before any symmetry copies are generated, so that every copy inherits
     * the same permutation and atom correspondence between copies is preserved.
     *
     * Only for managers whose inner loop accumulates into 1D distributions. The form-factor
     * managers write a Distribution3D of n_ff^2 x bin_count entries - 651 kB for 8 active
     * types and 434 bins - so they are far past the crossover, and their locality is
     * additionally helped by consecutive atoms sharing a form-factor type. Decorrelating
     * their order would spread the working set over the whole 3D distribution.
     *
     * @tparam weighted_bins Selects the histogram entry width, which sets the footprint.
     * @param bin_count The bin count the histograms will be allocated with.
     * @param sets Any number of coordinate sets, or (possibly nested) containers of them.
     */
    template<bool weighted_bins, typename... Sets>
    void decorrelate_order(unsigned int bin_count, Sets&... sets) {
        if (!atom_order::is_beneficial<weighted_bins>(bin_count)) {return;}
        unsigned int seed = 0x9e3779b9u; // fixed, so runs stay reproducible
        (atom_order::shuffle_all(sets, seed), ...);
    }
}
