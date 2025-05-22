#pragma once

#include <hist/distance_calculator/SimpleKernel.h>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief GPU kernel for simple histogram calculations based on WebGPU.
     */
    template<bool weighted_bins>
    struct WebGPUSimple : public SimpleKernel<weighted_bins> {
        WebGPUSimple();
        ~WebGPUSimple() override;

        void initialize();

        int enqueue_calculate_self(const hist::detail::CompactCoordinates& a, int scaling = 1, int merge_id = -1) override;
        int enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2, int scaling = 1, int merge_id = -1) override;
        SimpleKernel<weighted_bins>::run_result run() override;
    };
}