#pragma once

#include <gpu/WebGPU/WebGPU.h>
#include <gpu/WebGPU/GPUInstance.h>
#include <gpu/WebGPU/Buffers.h>
#include <hist/distance_calculator/SimpleKernel.h>
#include <constants/ConstantsAxes.h>
#include <utility/observer_ptr.h>

#include <unordered_map>

namespace ausaxs::gpu {
    template<bool weighted_bins>
    struct BufferManager {
        using GenericDistribution1D_t = typename ausaxs::hist::GenericDistribution1D<weighted_bins>::type;
        using run_result = typename ausaxs::hist::distance_calculator::SimpleKernel<weighted_bins>::run_result;

        int manage_self(wgpu::Buffer buffer, int merge_id = -1);
        int manage_cross(wgpu::Buffer buffer, int merge_id = -1);
        run_result merge(ausaxs::gpu::GPUInstance& instance) const;

        // vector of vectors to allow for multiple buffers per merge_id
        // this is different from the CPU version where only a single buffer per thread is needed, as they can simply accumulate results in a single vector
        std::vector<std::vector<wgpu::Buffer>> self_results, cross_results;
        std::unordered_map<int, int> self_merge_ids, cross_merge_ids;
    };
}