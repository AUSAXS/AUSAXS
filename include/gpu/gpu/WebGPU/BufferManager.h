#pragma once

#include <webgpu/webgpu.hpp>

#include <gpu/WebGPU/InstanceManager.h>
#include <constants/ConstantsAxes.h>
#include <utility/observer_ptr.h>

#include <unordered_map>

namespace ausaxs::gpu {
    struct BufferManager {
        struct run_result {
            std::unordered_map<int, std::vector<double>> self;
            std::unordered_map<int, std::vector<double>> cross;
        };

        int manage_self(wgpu::Buffer buffer, int merge_id = -1);
        int manage_cross(wgpu::Buffer buffer, int merge_id = -1);
        run_result merge(ausaxs::gpu::InstanceManager& instance) const;

        // vector of vectors to allow for multiple buffers per merge_id
        // this is different from the CPU version where only a single buffer per thread is needed, as they can simply accumulate results in a single vector
        std::vector<std::vector<wgpu::Buffer>> self_results, cross_results;
        std::unordered_map<int, int> self_merge_ids, cross_merge_ids;
    };
}

inline int ausaxs::gpu::BufferManager::manage_self(wgpu::Buffer buffer, int merge_id) {
    int res_idx;
    if (!self_merge_ids.contains(merge_id)) {
        res_idx = static_cast<int>(self_results.size());
        merge_id = merge_id == -1 ? res_idx : merge_id;
        self_merge_ids[merge_id] = res_idx;
        self_results.emplace_back(std::vector{buffer}); // new vector for the new merge_id
    } else {
        res_idx = self_merge_ids[merge_id];
        self_results[res_idx].emplace_back(buffer); // append to existing vector
    }
    return res_idx;
}

inline int ausaxs::gpu::BufferManager::manage_cross(wgpu::Buffer buffer, int merge_id) {
    int res_idx;
    if (!cross_merge_ids.contains(merge_id)) {
        res_idx = static_cast<int>(cross_results.size());
        merge_id = merge_id == -1 ? res_idx : merge_id;
        cross_merge_ids[merge_id] = res_idx;
        cross_results.emplace_back(std::vector{buffer}); // new vector for the new merge_id
    } else {
        res_idx = cross_merge_ids[merge_id];
        cross_results[res_idx].emplace_back(buffer); // append to existing vector
    }
    return res_idx;
}

inline void merge_buffer(ausaxs::gpu::InstanceManager& instance, std::vector<double>& destination, wgpu::Buffer histogram_readback) {
    assert(histogram_readback && "Readback buffer is null.");
    assert(histogram_readback.getUsage() & (wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead) && "Readback buffer does not have the correct usage flags.");

    struct CallbackInfo {
        ausaxs::observer_ptr<std::vector<double>> cpu;
        ausaxs::observer_ptr<wgpu::Buffer> gpu;
    } info;
    info.cpu = &destination;
    info.gpu = &histogram_readback;

    bool done = false;
    wgpu::BufferMapCallbackInfo map_callback;
    map_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
    map_callback.userdata1 = &info;
    map_callback.userdata2 = &done;
    map_callback.callback = [] (WGPUMapAsyncStatus, WGPUStringView, void* info_p, void* done_p) -> void {
        CallbackInfo& info = *reinterpret_cast<CallbackInfo*>(info_p);
        bool& done = *reinterpret_cast<bool*>(done_p);
        const float* output = (const float*) info.gpu->getConstMappedRange(0, info.gpu->getSize());
        std::transform(info.cpu->begin(), info.cpu->end(), output, info.cpu->begin(), std::plus<double>());
        info.gpu->unmap();
        done = true;
    };
    histogram_readback.mapAsync(wgpu::MapMode::Read, 0, histogram_readback.getSize(), map_callback);
    instance.process();
    assert(done && "Readback buffer mapping did not finish.");
}

inline std::vector<std::vector<wgpu::Buffer>> create_readback_buffers(const ausaxs::gpu::InstanceManager& instance, std::vector<std::vector<wgpu::Buffer>> hists) {
    wgpu::CommandEncoder encoder = instance.device.createCommandEncoder();

    static wgpu::BufferDescriptor readback_buffer_desc = [] () {
        wgpu::BufferDescriptor desc;
        desc.size = ausaxs::constants::axes::d_vals.size()*sizeof(float);
        desc.mappedAtCreation = false;
        desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
        return desc;
    }();

    std::vector<std::vector<wgpu::Buffer>> readback_buffers;
    for (const auto& vec : hists) {
        std::vector<wgpu::Buffer> buffers_i;
        for (const auto& buffer : vec) {
            assert(buffer && "Buffer is null.");
            wgpu::Buffer readback_buffer = instance.device.createBuffer(readback_buffer_desc);
            assert(readback_buffer.getSize() == buffer.getSize() && "Readback buffer size does not match histogram buffer size.");
            encoder.copyBufferToBuffer(buffer, 0, readback_buffer, 0, readback_buffer.getSize());
            buffers_i.emplace_back(readback_buffer);
        }
        readback_buffers.emplace_back(std::move(buffers_i));
    }

    auto command_buffer = encoder.finish();
    instance.device.getQueue().submit(command_buffer);
    command_buffer.release();
    encoder.release();
    return readback_buffers;
}

inline ausaxs::gpu::BufferManager::run_result ausaxs::gpu::BufferManager::merge(ausaxs::gpu::InstanceManager& instance) const {
    run_result result;

    auto self_readback_buffers = create_readback_buffers(instance, self_results);
    for (const auto& [merge_id, idx] : self_merge_ids) {
        std::vector<double> merged_hist(ausaxs::constants::axes::d_vals.size());
        auto& buffers = self_readback_buffers[idx];
        for (const auto& buffer : buffers) {
            merge_buffer(instance, merged_hist, buffer);
        }
        result.self[merge_id] = std::move(merged_hist);
    }

    auto cross_readback_buffers = create_readback_buffers(instance, cross_results);
    for (const auto& [merge_id, idx] : cross_merge_ids) {
        std::vector<double> merged_hist(ausaxs::constants::axes::d_vals.size());
        auto& buffers = self_readback_buffers[idx];
        for (const auto& buffer : buffers) {
            merge_buffer(instance, merged_hist, buffer);
        }
        result.cross[merge_id] = std::move(merged_hist);
    }

    return result;
}