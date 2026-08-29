// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/WebGPU/WebGPU.h>
#include <gpu/WebGPU/GPUInstance.h>
#include <gpu/WebGPU/Buffers.h>
#include <hist/distance_calculator/SimpleKernel.h>

#include <unordered_map>
#include <vector>

namespace ausaxs::gpu {
    /**
     * @brief Owns the histogram buffers of a single calculation.
     *
     * Unlike the CPU kernel, which needs one partial histogram per thread, all dispatches sharing a
     * merge_id accumulate atomically into a single GPU buffer. The manager therefore only has to hand
     * out the buffer belonging to a merge_id, and read one buffer per merge_id back at the end.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class BufferManager {
        public:
            using GenericDistribution1D_t = typename ausaxs::hist::GenericDistribution1D<weighted_bins>::type;
            using run_result = typename ausaxs::hist::distance_calculator::SimpleKernel<weighted_bins, variable_bin_width>::run_result;

            /**
             * @brief Get the histogram buffer to accumulate a self-correlation into, creating it if necessary.
             *
             * @param merge_id The requested merge id, or -1 to allocate a new one.
             * @return The buffer and the index of the calculation in the result.
             */
            std::pair<wgpu::Buffer, int> get_self_buffer(int merge_id);
            std::pair<wgpu::Buffer, int> get_cross_buffer(int merge_id); //< @copydoc get_self_buffer

            /**
             * @brief Copy all histogram buffers back to the host and convert them to distributions.
             *        This must only be called after the dispatches writing them have been submitted.
             */
            run_result merge();

            /// @brief Release all owned buffers and reset the manager to its initial state.
            void clear();

        private:
            std::vector<wgpu::Buffer> self_results, cross_results;
            std::unordered_map<int, int> self_merge_ids, cross_merge_ids;

            std::pair<wgpu::Buffer, int> get_buffer(
                int merge_id,
                std::vector<wgpu::Buffer>& results,
                std::unordered_map<int, int>& merge_ids
            );
    };
}
