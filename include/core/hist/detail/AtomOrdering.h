// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <settings/InternalState.h>

#include <ranges>

namespace ausaxs::hist::detail {
    namespace atom_order {
        /**
         * @brief The largest histogram footprint for which decorrelating the atom order pays off.
         *        This is only limited by the L1 cache size, remembering that it may be shared between two SMT threads. 
         */
        constexpr int max_histogram_bytes = 24*1024;

        /**
         * @brief The bytes one histogram of @a bin_count bins occupies.
         *        Weighted bins track a running centre and count per bin, so their entries are wider.
         */
        template<bool weighted_bins>
        constexpr int histogram_bytes(int bin_count) {
            if constexpr (weighted_bins) {
                return bin_count*static_cast<int>(sizeof(hist::detail::WeightedEntry));
            } else {
                return bin_count*static_cast<int>(sizeof(constants::axes::d_type));
            }
        }

        /**
         * @brief Whether decorrelating the atom order is expected to be a net win.
         */
        template<bool weighted_bins>
        constexpr bool is_beneficial(int bin_count) {
            return histogram_bytes<weighted_bins>(bin_count) <= max_histogram_bytes;
        }

        // a coordinate set that can be permuted directly
        template<typename T>
        concept Shufflable = requires(T& t) {t.shuffle_order();};

        template<Shufflable T>
        void shuffle_all(T& set) {set.shuffle_order();}

        // a container of sets - possibly nested, as in the per-body symmetry data
        template<std::ranges::input_range Range> requires (!Shufflable<Range>)
        void shuffle_all(Range& sets) {
            for (auto& set : sets) {shuffle_all(set);}
        }
    }

    /**
     * @brief Permute the atom order of the given coordinate sets, if the histogram is small enough for that to speed up the accumulation.
     */
    template<bool weighted_bins, typename... Sets>
    void decorrelate_order(int bin_count, Sets&... sets) {
        if (!settings::internal_state::allow_decorrelate_atom_order) {return;}
        if (!atom_order::is_beneficial<weighted_bins>(bin_count)) {return;}
        (atom_order::shuffle_all(sets), ...);
    }
}
