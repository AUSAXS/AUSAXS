// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/WebGPU/WebGPU.h>
#include <gpu/WebGPU/GPUInstance.h>
#include <gpu/WebGPU/Buffers.h>
#include <gpu/WebGPU/BufferManager.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/SimpleKernel.h>

#include <vector>

namespace ausaxs::gpu {
    /**
     * @brief Evaluates pairwise distance histograms with WebGPU compute shaders.
     *
     * Each submit call uploads its coordinates and records a dispatch, but no work is sent to the GPU
     * until run() is called. All dispatches sharing a merge_id accumulate atomically into the same
     * histogram buffer, so only one buffer per merge_id has to be read back.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class WebGPUBackend {
        public:
            using run_result = typename hist::distance_calculator::SimpleKernel<weighted_bins, variable_bin_width>::run_result;

            WebGPUBackend();
            ~WebGPUBackend();

            int submit_self(const hist::detail::CompactCoordinates<variable_bin_width>& atoms, int scaling = 1, int merge_id = -1);
            int submit_cross(
                const hist::detail::CompactCoordinates<variable_bin_width>& atoms_1,
                const hist::detail::CompactCoordinates<variable_bin_width>& atoms_2,
                int scaling = 1,
                int merge_id = -1
            );

            /// @brief Dispatch all submitted work, wait for it to finish and read back the histograms.
            run_result run();

            int size_self_result() const;
            int size_cross_result() const;

        private:
            BufferManager<weighted_bins, variable_bin_width> buffer_manager;
            std::vector<wgpu::Buffer> temporary_buffers;        // per-dispatch buffers, released after run()
            std::vector<wgpu::BindGroup> temporary_bind_groups; // per-dispatch bind groups, released after run()
            wgpu::CommandEncoder encoder = nullptr;             // accumulates the dispatches until run()
            int self_result_count = 0, cross_result_count = 0;

            wgpu::BindGroup assign_buffers(wgpu::Buffer atoms_1, wgpu::Buffer atoms_2, wgpu::Buffer histogram, wgpu::Buffer parameters);
            void dispatch(wgpu::BindGroup bind_group, wgpu::ComputePipeline pipeline, std::size_t atom_count);
            void release_temporaries();
    };
}
