// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/distance_calculator/detail/CalculatorCPU.h>
#include <hist/distance_calculator/detail/GPUKernel.h>
#include <settings/GeneralSettings.h>

#include <cassert>
#include <optional>
#include <span>
#include <vector>

namespace ausaxs::hist {
    constexpr bool TRACK_FF = true;
    constexpr bool WITHOUT_FF = false;
}

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Queues pairwise distance histogram calculations on the CPU or GPU kernel, into the results of a HistogramStore.
     *        The backend is picked once, on construction, from settings::general::gpu.
     *
     * A coordinate set is either flat, or partitioned into the store's classes() as one set per class. Each calculation
     * writes into a result of the store of the shape its sets call for, see HistogramStore. Several calculations may
     * write into the same result, in which case they are summed.
     *
     * Counting convention: a cross-correlation counts every pair of points pair_factor times. A self-correlation counts
     * every pair 2*scaling times, and every point with itself scaling times. The self-correlation of a partitioned set
     * puts each pair of unlike classes in the (k1, k2) histogram with k1 < k2 only, while the cross-correlation of two
     * partitioned sets fills both orderings, so a result mixing the two must be read symmetrically in the class pair.
     *
     * Each pair contributes the product of the weights of its two points, unless @a unit_weights is set: then every point
     * weighs 1, the stored weights are never read, and the histograms are plain pair counts. This is for calculations whose
     * weighting is applied later, by class, such as the form factors applied by the intensity calculator.
     *
     * All data references must be valid until run() is called.
     */
    template<bool weighted_bins, bool variable_bin_width, bool unit_weights>
    class Calculator {
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        using PartitionedCoordinates_t = std::vector<CompactCoordinates_t>;
        using Row = std::span<typename HistogramStore<weighted_bins>::entry_type>;
        public:
            /**
             * @brief Construct a calculator writing into @a store, which must outlive it.
             */
            explicit Calculator(HistogramStore<weighted_bins>& store) : store(&store) {
                if (settings::general::gpu) {gpu.emplace(store);}
                else {cpu.emplace(store);}
            }

            /**
             * @brief Queue the self-correlation of the flat set @a a into the allocate_1d() @a result.
             */
            void enqueue_calculate_self(const CompactCoordinates_t& a, int result, int scaling = 1) {
                self(a, store->row(result), scaling);
            }

            /**
             * @brief Queue the self-correlation of the partitioned set @a a into the allocate_3d() @a result.
             */
            void enqueue_calculate_self(const PartitionedCoordinates_t& a, int result, int scaling = 1) {
                assert(static_cast<int>(a.size()) == store->classes() && "Calculator: expected one set per class of the store.");
                for (int k1 = 0; k1 < static_cast<int>(a.size()); ++k1) {
                    self(a[k1], store->row(result, k1, k1), scaling);
                    for (int k2 = k1+1; k2 < static_cast<int>(a.size()); ++k2) {
                        cross(a[k1], a[k2], store->row(result, k1, k2), 2*scaling);
                    }
                }
            }

            /**
             * @brief Queue the cross-correlation of the flat sets @a a1 and @a a2 into the allocate_1d() @a result.
             */
            void enqueue_calculate_cross(const CompactCoordinates_t& a1, const CompactCoordinates_t& a2, int result, int pair_factor) {
                cross(a1, a2, store->row(result), pair_factor);
            }

            /**
             * @brief Queue the cross-correlation of the partitioned set @a a1 and the flat set @a a2 into the allocate_2d() @a result.
             */
            void enqueue_calculate_cross(const PartitionedCoordinates_t& a1, const CompactCoordinates_t& a2, int result, int pair_factor) {
                assert(static_cast<int>(a1.size()) == store->classes() && "Calculator: expected one set per class of the store.");
                for (int k = 0; k < static_cast<int>(a1.size()); ++k) {
                    cross(a1[k], a2, store->row(result, k), pair_factor);
                }
            }

            /**
             * @brief Queue the cross-correlation of the partitioned sets @a a1 and @a a2 into the allocate_3d() @a result.
             *        The (k1, k2) histogram holds class k1 of @a a1 against class k2 of @a a2.
             */
            void enqueue_calculate_cross(const PartitionedCoordinates_t& a1, const PartitionedCoordinates_t& a2, int result, int pair_factor) {
                assert(static_cast<int>(a1.size()) == store->classes() && static_cast<int>(a2.size()) == store->classes() && "Calculator: expected one set per class of the store.");
                for (int k1 = 0; k1 < static_cast<int>(a1.size()); ++k1) {
                    for (int k2 = 0; k2 < static_cast<int>(a2.size()); ++k2) {
                        cross(a1[k1], a2[k2], store->row(result, k1, k2), pair_factor);
                    }
                }
            }

            /**
             * @brief Withhold everything enqueued from here on, as one group, until release_hold().
             *
             * This only does anything on the GPU backend, where each submission would otherwise be dispatched to the device
             * on its own; holding lets a group of related jobs go as one. The CPU kernel dispatches every job immediately regardless.
             */
            void hold() {if (gpu) {gpu->hold();}}

            /**
             * @brief Dispatch the groups held since hold(), and stop holding.
             */
            void release_hold() {if (gpu) {gpu->release_hold();}}

            /**
             * @brief Calculate the queued histograms into their results in the store.
             *        This will block until all calculations are done.
             */
            void run() {
                assert((cpu.has_value() || gpu.has_value()) && "Calculator: the constructor engages exactly one backend.");
                if (cpu) {cpu->run();}
                else {gpu->run();}
            }

        private:
            observer_ptr<HistogramStore<weighted_bins>> store;

            // exactly one of these is engaged, as decided by the constructor
            std::optional<detail::CalculatorCPU<weighted_bins, variable_bin_width, unit_weights>> cpu;
            std::optional<detail::GPUKernel<weighted_bins, variable_bin_width, unit_weights>> gpu;

            void self(const CompactCoordinates_t& a, Row row, int scaling) {
                assert((cpu.has_value() || gpu.has_value()) && "Calculator: the constructor engages exactly one backend.");
                if (cpu) {cpu->enqueue_calculate_self(a, row, scaling);}
                else {gpu->enqueue_calculate_self(a, row, scaling);}
            }

            void cross(const CompactCoordinates_t& a1, const CompactCoordinates_t& a2, Row row, int pair_factor) {
                assert((cpu.has_value() || gpu.has_value()) && "Calculator: the constructor engages exactly one backend.");
                if (cpu) {cpu->enqueue_calculate_cross(a1, a2, row, pair_factor);}
                else {gpu->enqueue_calculate_cross(a1, a2, row, pair_factor);}
            }
    };
}
