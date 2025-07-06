#pragma once

#include <hist/distance_calculator/SimpleKernel.h>
#include <gpu/WebGPU/WebGPUController.h>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief GPU kernel for simple histogram calculations based on WebGPU.
     */
    template<bool weighted_bins>
    class WebGPUSimple : public SimpleKernel<weighted_bins> {
        public:
            WebGPUSimple();
            ~WebGPUSimple() override;

            int enqueue_calculate_self(const hist::detail::CompactCoordinates& a, int scaling = 1, int merge_id = -1) override;
            int enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2, int scaling = 1, int merge_id = -1) override;
            SimpleKernel<weighted_bins>::run_result run() override;
        
        private:
            gpu::WebGPU<weighted_bins> controller;
    };
}

template<bool weighted_bins>
inline int ausaxs::hist::distance_calculator::WebGPUSimple<weighted_bins>::enqueue_calculate_self(const hist::detail::CompactCoordinates& a, int scaling, int merge_id) {
    return controller.submit_self(a, merge_id);
}

template<bool weighted_bins>
inline int ausaxs::hist::distance_calculator::WebGPUSimple<weighted_bins>::enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2, int scaling, int merge_id) {
    return controller.submit_cross(a1, a2, merge_id);
}

template<bool weighted_bins>
inline typename ausaxs::hist::distance_calculator::SimpleKernel<weighted_bins>::run_result ausaxs::hist::distance_calculator::WebGPUSimple<weighted_bins>::run() {
    return controller.run();
}