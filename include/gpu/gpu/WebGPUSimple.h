// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distance_calculator/SimpleKernel.h>
#include <gpu/WebGPU/WebGPUBackend.h>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief WebGPU kernel for simple histogram calculations. See gpu::WebGPUBackend for the implementation.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class WebGPUSimple : public SimpleKernel<weighted_bins, variable_bin_width> {
        public:
            using run_result = typename SimpleKernel<weighted_bins, variable_bin_width>::run_result;

            int enqueue_calculate_self(
                const hist::detail::CompactCoordinates<variable_bin_width>& a,
                int scaling = 1,
                int merge_id = -1
            ) override {
                return backend.submit_self(a, scaling, merge_id);
            }

            int enqueue_calculate_cross(
                const hist::detail::CompactCoordinates<variable_bin_width>& a1,
                const hist::detail::CompactCoordinates<variable_bin_width>& a2,
                int scaling = 1,
                int merge_id = -1
            ) override {
                return backend.submit_cross(a1, a2, scaling, merge_id);
            }

            int size_self_result() const override {return backend.size_self_result();}
            int size_cross_result() const override {return backend.size_cross_result();}

            run_result run() override {return backend.run();}

        private:
            gpu::WebGPUBackend<weighted_bins, variable_bin_width> backend;
    };
}
