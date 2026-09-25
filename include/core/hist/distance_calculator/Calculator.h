// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/HistogramStore.h>
#include <hist/distance_calculator/detail/CalculatorCPU.h>
#include <hist/distance_calculator/detail/GPUKernel.h>
#include <settings/GeneralSettings.h>

#include <cassert>
#include <optional>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Queues pairwise distance histogram calculations on the CPU or GPU kernel, into the rows of a HistogramStore.
     *        The backend is picked once, on construction, from settings::general::gpu.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class Calculator {
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        public:
            /**
             * @brief Construct a calculator accumulating into @a store, which must outlive it.
             */
            explicit Calculator(HistogramStore<weighted_bins>& store) {
                if (settings::general::gpu) {gpu.emplace(store);}
                else {cpu.emplace(store);}
            }

            /**
             * @brief Queue a self-correlation calculation.
             *        The data reference must be valid until run() is called. 
             *
             * @param a The data to calculate the self-correlation for. The reference must be valid until run() is called.
             * @param h The handle to accumulate into.
             * @param scaling The scaling factor to apply to the pair counts. 
             */
            void enqueue_calculate_self(const CompactCoordinates_t& a, int h, int scaling = 1) {
                assert((cpu.has_value() || gpu.has_value()) && "Calculator: the constructor engages exactly one backend.");
                if (cpu) {cpu->enqueue_calculate_self(a, h, scaling);}
                else {gpu->enqueue_calculate_self(a, h, scaling);}
            }

            /**
             * @brief Queue a cross-correlation calculation.
             *        The data references must be valid until run() is called. 
             *
             * @param a1 The first set of data to calculate the cross-correlation for.
             * @param a2 The second set of data to calculate the cross-correlation for.
             * @param h The handle to accumulate into.
             * @param scaling The scaling factor to apply to the pair counts. 
             */
            void enqueue_calculate_cross(const CompactCoordinates_t& a1, const CompactCoordinates_t& a2, int h, int scaling) {
                assert((cpu.has_value() || gpu.has_value()) && "Calculator: the constructor engages exactly one backend.");
                if (cpu) {cpu->enqueue_calculate_cross(a1, a2, h, scaling);}
                else {gpu->enqueue_calculate_cross(a1, a2, h, scaling);}
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
             * @brief Calculate the queued histograms into their rows of the store.
             *        This will block until all calculations are done.
             */
            void run() {
                assert((cpu.has_value() || gpu.has_value()) && "Calculator: the constructor engages exactly one backend.");
                if (cpu) {cpu->run();}
                else {gpu->run();}
            }

        private:
            // exactly one of these is engaged, as decided by the constructor
            std::optional<detail::CalculatorCPU<weighted_bins, variable_bin_width>> cpu;
            std::optional<detail::GPUKernel<weighted_bins, variable_bin_width>> gpu;
    };
}
