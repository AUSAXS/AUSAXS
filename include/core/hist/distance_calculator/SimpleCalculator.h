// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/SimpleCPU.h>
#include <hist/distance_calculator/SimpleGPU.h>
#include <settings/GeneralSettings.h>

#include <cassert>
#include <optional>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Queues pairwise distance histogram calculations on the CPU or GPU kernel.
     *        The backend is picked once, on construction, from settings::general::gpu.
     *
     * The caller must keep all submitted data alive until run() returns.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class SimpleCalculator {
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        public:
            using run_result = typename SimpleCPU<weighted_bins, variable_bin_width>::run_result;

            /**
             * @brief Construct a calculator whose result histograms span @a bin_count bins.
             */
            explicit SimpleCalculator(int bin_count) {
                if (settings::general::gpu) {gpu.emplace(bin_count);}
                else {cpu.emplace(bin_count);}
            }

            /**
             * @brief Queue a self-correlation calculation.
             *        This is faster than calling the cross-correlation method with the same data, as some optimizations can be made.
             *
             * @param a The data to calculate the self-correlation for. The reference must be valid until run() is called.
             * @param scaling The scaling factor to apply to the result.
             * @param merge_id The result vector id this calculation can be merged into. Supplying this can save significant memory resources.
             *
             * @return The index of the data in the result vector.
             */
            int enqueue_calculate_self(const CompactCoordinates_t& a, int scaling = 1, int merge_id = -1) {
                assert((cpu.has_value() || gpu.has_value()) && "SimpleCalculator: the constructor engages exactly one backend.");
                return cpu ? cpu->enqueue_calculate_self(a, scaling, merge_id)
                           : gpu->enqueue_calculate_self(a, scaling, merge_id);
            }

            /**
             * @brief Queue a cross-correlation calculation.
             *
             * @param a1 The first set of data to calculate the cross-correlation for. The reference must be valid until run() is called.
             * @param a2 The second set of data to calculate the cross-correlation for. The reference must be valid until run() is called.
             * @param scaling The scaling factor to apply to the result.
             * @param merge_id The result vector id this calculation can be merged into. Supplying this can save significant memory resources.
             *
             * @return The index of the data in the result vector.
             */
            int enqueue_calculate_cross(const CompactCoordinates_t& a1, const CompactCoordinates_t& a2, int scaling = 1, int merge_id = -1) {
                assert((cpu.has_value() || gpu.has_value()) && "SimpleCalculator: the constructor engages exactly one backend.");
                return cpu ? cpu->enqueue_calculate_cross(a1, a2, scaling, merge_id)
                           : gpu->enqueue_calculate_cross(a1, a2, scaling, merge_id);
            }

            /**
             * @brief Withhold everything enqueued from here on, as one group, until release_hold().
             *
             * This is useful for batching calculations that share a merge_id, so they can be dispatched together. This only does
             * anything on the GPU backend, where each submission would otherwise be allocated its own buffer on the device;
             * holding allows them to share. The CPU kernel dispatches every job immediately regardless.
             */
            void hold() {if (gpu) {gpu->hold();}}

            /**
             * @brief Dispatch the group held since hold(), and stop holding.
             *        Call this as soon as the group is complete, as the work proceeds asynchronously from
             *        there: whatever the caller does before run() overlaps with it.
             */
            void release_hold() {if (gpu) {gpu->release_hold();}}

            /**
             * @brief Get the current size of the result vector.
             */
            int size_self_result() const {assert((cpu.has_value() || gpu.has_value()) && "SimpleCalculator: the constructor engages exactly one backend."); return cpu ? cpu->size_self_result() : gpu->size_self_result();}
            int size_cross_result() const {assert((cpu.has_value() || gpu.has_value()) && "SimpleCalculator: the constructor engages exactly one backend."); return cpu ? cpu->size_cross_result() : gpu->size_cross_result();} //< @copydoc size_self_result

            /**
             * @brief Calculate the queued histograms.
             *        This will block until all calculations are done.
             *
             * @return The calculated histograms.
             */
            run_result run() {assert((cpu.has_value() || gpu.has_value()) && "SimpleCalculator: the constructor engages exactly one backend."); return cpu ? cpu->run() : gpu->run();}

        private:
            // exactly one of these is engaged, as decided by the constructor
            std::optional<SimpleCPU<weighted_bins, variable_bin_width>> cpu;
            std::optional<SimpleGPU<weighted_bins, variable_bin_width>> gpu;
    };
}
