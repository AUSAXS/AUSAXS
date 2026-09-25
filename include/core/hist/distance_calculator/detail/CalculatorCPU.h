// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/distance_calculator/detail/CPUKernel.h>
#include <utility/MultiThreading.h>
#include <utility/observer_ptr.h>

namespace ausaxs::hist::distance_calculator::detail {
    /**
     * @brief Queues and evaluates pairwise distance histograms on the global thread pool, into the rows of a HistogramStore.
     *        The caller must keep all submitted data alive until run() returns.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class CalculatorCPU {
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        public:
            /**
             * @brief Construct a calculator accumulating into @a store, which must outlive it.
             */
            explicit CalculatorCPU(HistogramStore<weighted_bins>& store) : store(&store) {}

            /**
             * @brief Queue the self-correlation of @a a into the row @a h: each pair counted 2*@a scaling times, and each
             *        point with itself @a scaling times.
             */
            void enqueue_calculate_self(const CompactCoordinates_t& a, int h, int scaling) {
                if (a.empty()) {store->clear(h); return;}
                dispatch_scaling(2*scaling, [&a, target = store->target(h)] (auto pair) {
                    constexpr int pair_factor = decltype(pair)::value;
                    if constexpr (pair_factor % 2 == 0) { // always the case, since the factor is twice the scaling
                        enqueue_self<pair_factor, pair_factor/2>(a, target);
                    }
                });
            }

            /**
             * @brief Queue the cross-correlation of @a a1 and @a a2 into the row @a h, each pair counted @a pair_factor times.
             *        The work is chunked over whichever of the two is larger.
             */
            void enqueue_calculate_cross(const CompactCoordinates_t& a1, const CompactCoordinates_t& a2, int h, int pair_factor) {
                if (a1.empty() || a2.empty()) {store->clear(h); return;}
                dispatch_scaling(pair_factor, [&a1, &a2, target = store->target(h)] (auto pair) {
                    enqueue_balanced_cross<decltype(pair)::value>(a1, a2, target);
                });
            }

            /**
             * @brief Wait for the queued calculations. This will block until all calculations are done.
             */
            void run() {
                utility::multi_threading::get_global_pool()->wait();
                store->fold();
            }

        private:
            observer_ptr<HistogramStore<weighted_bins>> store;
    };
}
