// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/WebGPU/BufferManager.h>
#include <settings/HistogramSettings.h>
#include <utility/observer_ptr.h>

#include <algorithm>
#include <cassert>

using namespace ausaxs;
using namespace ausaxs::gpu;

template<bool weighted_bins, bool variable_bin_width>
std::pair<wgpu::Buffer, int> BufferManager<weighted_bins, variable_bin_width>::get_buffer(
    int merge_id,
    std::vector<wgpu::Buffer>& results,
    std::unordered_map<int, int>& merge_ids
) {
    if (auto it = merge_ids.find(merge_id); merge_id != -1 && it != merge_ids.end()) {
        return {results[it->second], it->second};
    }

    // no buffer is accumulating for this merge_id yet, so create one
    int res_idx = static_cast<int>(results.size());
    merge_ids[merge_id == -1 ? res_idx : merge_id] = res_idx;
    results.emplace_back(Buffers<weighted_bins>::create_histogram_buffer());
    return {results.back(), res_idx};
}

template<bool weighted_bins, bool variable_bin_width>
std::pair<wgpu::Buffer, int> BufferManager<weighted_bins, variable_bin_width>::get_self_buffer(int merge_id) {
    return get_buffer(merge_id, self_results, self_merge_ids);
}

template<bool weighted_bins, bool variable_bin_width>
std::pair<wgpu::Buffer, int> BufferManager<weighted_bins, variable_bin_width>::get_cross_buffer(int merge_id) {
    return get_buffer(merge_id, cross_results, cross_merge_ids);
}

namespace {
    // copy the histograms into mappable buffers. this is done in a single command buffer for all of them.
    template<bool weighted_bins>
    std::vector<wgpu::Buffer> create_readback_buffers(const std::vector<wgpu::Buffer>& histograms) {
        if (histograms.empty()) {return {};}
        auto& gpu_instance = GPUInstance::get();
        wgpu::CommandEncoder encoder = gpu_instance.device.createCommandEncoder();

        std::vector<wgpu::Buffer> readback_buffers;
        readback_buffers.reserve(histograms.size());
        for (const auto& buffer : histograms) {
            assert(buffer && "BufferManager: histogram buffer is null.");
            wgpu::Buffer readback_buffer = Buffers<weighted_bins>::create_readback_buffer();
            assert(readback_buffer.getSize() == buffer.getSize() && "BufferManager: readback buffer size does not match histogram buffer size.");
            encoder.copyBufferToBuffer(buffer, 0, readback_buffer, 0, readback_buffer.getSize());
            readback_buffers.emplace_back(readback_buffer);
        }

        auto command_buffer = encoder.finish();
        gpu_instance.queue.submit(command_buffer);
        command_buffer.release();
        encoder.release();
        return readback_buffers;
    }

    // map a readback buffer and convert its contents to the host distribution type
    template<bool weighted_bins>
    typename ausaxs::hist::GenericDistribution1D<weighted_bins>::type read_buffer(wgpu::Buffer readback) {
        using GPUHistogramType = typename Buffers<weighted_bins>::HistogramType;
        assert(readback && "BufferManager: readback buffer is null.");

        std::vector<GPUHistogramType> gpu_histogram(settings::axes::bin_count);
        assert(readback.getSize() == gpu_histogram.size()*sizeof(GPUHistogramType) && "BufferManager: readback buffer size does not match histogram size.");

        {   // read the GPU histogram into the CPU histogram
            struct CallbackInfo {
                observer_ptr<std::vector<GPUHistogramType>> cpu;
                observer_ptr<wgpu::Buffer> gpu;
            } info;
            info.cpu = &gpu_histogram;
            info.gpu = &readback;

            bool done = false;
            wgpu::BufferMapCallbackInfo map_callback;
            map_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
            map_callback.userdata1 = &info;
            map_callback.userdata2 = &done;
            map_callback.callback = [] (WGPUMapAsyncStatus status, WGPUStringView, void* info_p, void* done_p) -> void {
                CallbackInfo& info = *reinterpret_cast<CallbackInfo*>(info_p);
                bool& done = *reinterpret_cast<bool*>(done_p);
                if (status == static_cast<WGPUMapAsyncStatus>(wgpu::MapAsyncStatus::Success)) {
                    const GPUHistogramType* output = static_cast<const GPUHistogramType*>(info.gpu->getConstMappedRange(0, info.gpu->getSize()));
                    std::copy(output, output + info.cpu->size(), info.cpu->begin());
                    info.gpu->unmap();
                }
                done = true;
            };
            readback.mapAsync(wgpu::MapMode::Read, 0, readback.getSize(), map_callback);
            GPUInstance::get().wait(done);
        }

        typename ausaxs::hist::GenericDistribution1D<weighted_bins>::type histogram(settings::axes::bin_count);
        for (int i = 0; i < static_cast<int>(gpu_histogram.size()); ++i) {
            const auto& entry = gpu_histogram[i];
            if constexpr (weighted_bins) {
                histogram.add_index(i, ausaxs::hist::detail::WeightedEntry{entry.value, entry.count, entry.bin_center});
            } else {
                histogram.add_index(i, entry.value);
            }
        }
        return histogram;
    }
}

template<bool weighted_bins, bool variable_bin_width>
typename BufferManager<weighted_bins, variable_bin_width>::run_result BufferManager<weighted_bins, variable_bin_width>::merge() {
    run_result result;

    auto self_readback_buffers = create_readback_buffers<weighted_bins>(self_results);
    auto cross_readback_buffers = create_readback_buffers<weighted_bins>(cross_results);

    for (const auto& [merge_id, idx] : self_merge_ids) {
        result.self[merge_id] = read_buffer<weighted_bins>(self_readback_buffers[idx]);
    }

    for (const auto& [merge_id, idx] : cross_merge_ids) {
        result.cross[merge_id] = read_buffer<weighted_bins>(cross_readback_buffers[idx]);
    }

    std::for_each(self_readback_buffers.begin(), self_readback_buffers.end(), [] (wgpu::Buffer& b) {b.release();});
    std::for_each(cross_readback_buffers.begin(), cross_readback_buffers.end(), [] (wgpu::Buffer& b) {b.release();});
    return result;
}

template<bool weighted_bins, bool variable_bin_width>
void BufferManager<weighted_bins, variable_bin_width>::clear() {
    std::for_each(self_results.begin(), self_results.end(), [] (wgpu::Buffer& b) {b.destroy(); b.release();});
    std::for_each(cross_results.begin(), cross_results.end(), [] (wgpu::Buffer& b) {b.destroy(); b.release();});
    self_results.clear();
    cross_results.clear();
    self_merge_ids.clear();
    cross_merge_ids.clear();
}

template class ausaxs::gpu::BufferManager<false, false>;
template class ausaxs::gpu::BufferManager<false, true>;
template class ausaxs::gpu::BufferManager<true, false>;
template class ausaxs::gpu::BufferManager<true, true>;
