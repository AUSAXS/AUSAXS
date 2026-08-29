// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinates.h>

#include <unordered_map>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Interface for queueing pairwise distance histogram calculations.
     *
     * Data is submitted through the enqueue methods and evaluated no later than the call to run().
     * Whether the work is dispatched immediately or deferred until run() is up to the implementation;
     * the only guarantee is that all submitted data must be kept alive by the caller until run() returns.
     *
     * Jobs sharing a @c merge_id accumulate into the same result histogram, which saves memory when many
     * calculations contribute to a single histogram (e.g. symmetry copies). The integer @c scaling factor
     * multiplies a job's contribution; only a bounded set of values is supported (see the implementations).
     *
     * See SimpleCPU for the thread pool implementation, and SimpleCalculator for how one is chosen.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class SimpleKernel {
        using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
        public:
            struct run_result {
                std::unordered_map<int, GenericDistribution1D_t> self;
                std::unordered_map<int, GenericDistribution1D_t> cross;
            };

            virtual ~SimpleKernel() = default;

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
            virtual int enqueue_calculate_self(const hist::detail::CompactCoordinates<variable_bin_width>& a, int scaling = 1, int merge_id = -1) = 0;

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
            virtual int enqueue_calculate_cross(
                const hist::detail::CompactCoordinates<variable_bin_width>& a1,
                const hist::detail::CompactCoordinates<variable_bin_width>& a2,
                int scaling = 1,
                int merge_id = -1
            ) = 0;

            /**
             * @brief Get the current size of the result vector.
             */
            virtual int size_self_result() const = 0;
            virtual int size_cross_result() const = 0; //< @copydoc size_self_result

            /**
             * @brief Calculate the queued histograms.
             *        This will block until all calculations are done.
             *
             * @return The calculated histograms.
             */
            virtual run_result run() = 0;
    };
}
