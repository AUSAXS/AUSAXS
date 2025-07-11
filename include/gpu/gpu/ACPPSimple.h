#pragma once

#include <hist/distance_calculator/SimpleKernel.h>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief GPU kernel for simple histogram calculations based on AdaptiveCPP.
     */
    template<bool weighted_bins>
    class ACPPSimple : public SimpleKernel<weighted_bins> {
        int enqueue_calculate_self(const hist::detail::CompactCoordinates& a, int scaling = 1, int merge_id = -1) override;
        int enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2, int scaling = 1, int merge_id = -1) override;
        run_result run() override;
    };
}

template<bool weighted_bins>
int ausaxs::hist::distance_calculator::ACPPSimple<weighted_bins>::enqueue_calculate_self(const hist::detail::CompactCoordinates& a, int scaling = 1, int merge_id = -1) {
    return 0;
}

template<bool weighted_bins>
int ausaxs::hist::distance_calculator::ACPPSimple<weighted_bins>::enqueue_calculate_cross(const hist::detail::CompactCoordinates& a1, const hist::detail::CompactCoordinates& a2, int scaling = 1, int merge_id = -1) {
    return 0;
}

template<bool weighted_bins>
ausaxs::hist::distance_calculator::ACPPSimple<weighted_bins>::run_result ausaxs::hist::distance_calculator::ACPPSimple<weighted_bins>::run() {
    return {};
}