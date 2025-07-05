#pragma once

#include <webgpu/webgpu.hpp>

#include <unordered_map>

namespace ausaxs::gpu {
    struct BufferManager {
        struct run_result {
            std::unordered_map<int, std::vector<double>> self;
            std::unordered_map<int, std::vector<double>> cross;
        };

        int manage_self(wgpu::Buffer buffer, int merge_id = -1);
        int manage_cross(wgpu::Buffer buffer, int merge_id = -1);
        run_result merge() const;

        std::vector<wgpu::Buffer> self_results, cross_results;
        std::unordered_map<int, int> self_merge_ids, cross_merge_ids;
    };
}

inline int ausaxs::gpu::BufferManager::manage_self(wgpu::Buffer buffer, int merge_id) {
    int res_idx;
    if (!self_merge_ids.contains(merge_id)) {
        res_idx = static_cast<int>(self_results.size());
        merge_id = merge_id == -1 ? res_idx : merge_id;
        self_merge_ids[merge_id] = res_idx;
        self_results.emplace_back(buffer);
    } else {
        res_idx = self_merge_ids[merge_id];
        assert(self_results[res_idx].getSize() == buffer.getSize() && "The result buffer has the wrong size.");
    }
    return res_idx;
}

inline int ausaxs::gpu::BufferManager::manage_cross(wgpu::Buffer buffer, int merge_id) {
    int res_idx;
    if (!cross_merge_ids.contains(merge_id)) {
        res_idx = static_cast<int>(cross_results.size());
        merge_id = merge_id == -1 ? res_idx : merge_id;
        cross_merge_ids[merge_id] = res_idx;
        cross_results.emplace_back(buffer);
    } else {
        res_idx = cross_merge_ids[merge_id];
        assert(cross_results[res_idx].getSize() == buffer.getSize() && "The result buffer has the wrong size.");
    }
    return res_idx;
}