// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/WebGPU/WebGPUBackend.h>
#include <gpu/WebGPU/shaders/ShaderStorage.h>
#include <hist/detail/data/WidthControllers.h>
#include <settings/HistogramSettings.h>
#include <utility/Logging.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

using namespace ausaxs;
using namespace ausaxs::gpu;

namespace {
    constexpr std::size_t workgroup_size = 256; // must match @workgroup_size in the wgsl shaders

    std::size_t workgroups_for(std::size_t invocations) {
        return (invocations + workgroup_size - 1)/workgroup_size;
    }

    template<bool variable_bin_width>
    ShaderParameters make_parameters(int scaling) {
        return {
            .inv_width = hist::detail::WidthController<variable_bin_width>::get_inv_width(),
            .scale = static_cast<std::uint32_t>(scaling),
            .bin_count = settings::axes::bin_count
        };
    }
}

template<bool weighted_bins, bool variable_bin_width>
WebGPUBackend<weighted_bins, variable_bin_width>::WebGPUBackend() {
    // touch the device and shaders up front so any initialization cost is not attributed to the first dispatch
    GPUInstance::get();
    shader::Simple::get<weighted_bins>();
}

template<bool weighted_bins, bool variable_bin_width>
WebGPUBackend<weighted_bins, variable_bin_width>::~WebGPUBackend() {
    if (encoder) {encoder.release();}
    release_temporaries();
    buffer_manager.clear();
}

template<bool weighted_bins, bool variable_bin_width>
wgpu::BindGroup WebGPUBackend<weighted_bins, variable_bin_width>::assign_buffers(
    wgpu::Buffer atoms_1,
    wgpu::Buffer atoms_2,
    wgpu::Buffer histogram,
    wgpu::Buffer parameters
) {
    assert(atoms_1 && atoms_2 && histogram && parameters && "WebGPUBackend::assign_buffers: buffer is null.");

    std::array<wgpu::BindGroupEntry, 4> entries = {wgpu::Default, wgpu::Default, wgpu::Default, wgpu::Default};
    entries[0].binding = 0;
    entries[0].buffer = atoms_1;
    entries[0].size = atoms_1.getSize();

    entries[1].binding = 1;
    entries[1].buffer = atoms_2;
    entries[1].size = atoms_2.getSize();

    entries[2].binding = 2;
    entries[2].buffer = histogram;
    entries[2].size = histogram.getSize();

    entries[3].binding = 3;
    entries[3].buffer = parameters;
    entries[3].size = parameters.getSize();

    wgpu::BindGroupDescriptor bind_group_desc;
    bind_group_desc.layout = shader::Simple::get<weighted_bins>().bind_group_layout;
    bind_group_desc.entryCount = entries.size();
    bind_group_desc.entries = entries.data();
    return GPUInstance::get().device.createBindGroup(bind_group_desc);
}

template<bool weighted_bins, bool variable_bin_width>
void WebGPUBackend<weighted_bins, variable_bin_width>::dispatch(wgpu::BindGroup bind_group, wgpu::ComputePipeline pipeline, std::size_t workgroups) {
    if (workgroups == 0) {return;}
    if (!encoder) {encoder = GPUInstance::get().device.createCommandEncoder();}

    auto compute_pass = encoder.beginComputePass(wgpu::Default);
    compute_pass.setPipeline(pipeline);
    compute_pass.setBindGroup(0, bind_group, 0, nullptr);
    compute_pass.dispatchWorkgroups(workgroups, 1, 1);
    compute_pass.end();
    compute_pass.release();
}

template<bool weighted_bins, bool variable_bin_width>
int WebGPUBackend<weighted_bins, variable_bin_width>::submit_self(
    const hist::detail::CompactCoordinates<variable_bin_width>& atoms,
    int scaling,
    int merge_id
) {
    auto atom_buffer = Buffers<weighted_bins>::template create_atom_buffer<variable_bin_width>(atoms);
    auto parameter_buffer = Buffers<weighted_bins>::create_parameter_buffer(make_parameters<variable_bin_width>(scaling));
    auto [histogram, index] = buffer_manager.get_self_buffer(merge_id);
    self_result_count = std::max(self_result_count, index+1);

    // each invocation evaluates two rows of the distance matrix, one from either end
    auto bind_group = assign_buffers(atom_buffer, Buffers<weighted_bins>::dummy_buffer(), histogram, parameter_buffer);
    dispatch(bind_group, shader::Simple::get<weighted_bins>().pipelines.self, workgroups_for((atoms.size() + 1)/2));

    temporary_buffers.emplace_back(atom_buffer);
    temporary_buffers.emplace_back(parameter_buffer);
    temporary_bind_groups.emplace_back(bind_group);
    return index;
}

template<bool weighted_bins, bool variable_bin_width>
int WebGPUBackend<weighted_bins, variable_bin_width>::submit_cross(
    const hist::detail::CompactCoordinates<variable_bin_width>& atoms_1,
    const hist::detail::CompactCoordinates<variable_bin_width>& atoms_2,
    int scaling,
    int merge_id
) {
    auto atom_buffer_1 = Buffers<weighted_bins>::template create_atom_buffer<variable_bin_width>(atoms_1);
    auto atom_buffer_2 = Buffers<weighted_bins>::template create_atom_buffer<variable_bin_width>(atoms_2);
    auto parameter_buffer = Buffers<weighted_bins>::create_parameter_buffer(make_parameters<variable_bin_width>(scaling));
    auto [histogram, index] = buffer_manager.get_cross_buffer(merge_id);
    cross_result_count = std::max(cross_result_count, index+1);

    auto bind_group = assign_buffers(atom_buffer_1, atom_buffer_2, histogram, parameter_buffer);
    dispatch(bind_group, shader::Simple::get<weighted_bins>().pipelines.cross, workgroups_for(atoms_1.size()));

    temporary_buffers.emplace_back(atom_buffer_1);
    temporary_buffers.emplace_back(atom_buffer_2);
    temporary_buffers.emplace_back(parameter_buffer);
    temporary_bind_groups.emplace_back(bind_group);
    return index;
}

template<bool weighted_bins, bool variable_bin_width>
typename WebGPUBackend<weighted_bins, variable_bin_width>::run_result WebGPUBackend<weighted_bins, variable_bin_width>::run() {
    if (encoder) { // submit the recorded dispatches
        auto command_buffer = encoder.finish();
        GPUInstance::get().queue.submit(command_buffer);
        command_buffer.release();
        encoder.release();
        encoder = nullptr;
    }

    auto result = buffer_manager.merge(); // blocks until the dispatches above have completed

    release_temporaries();
    buffer_manager.clear();
    self_result_count = 0;
    cross_result_count = 0;
    return result;
}

template<bool weighted_bins, bool variable_bin_width>
void WebGPUBackend<weighted_bins, variable_bin_width>::release_temporaries() {
    std::for_each(temporary_bind_groups.begin(), temporary_bind_groups.end(), [] (wgpu::BindGroup& b) {b.release();});
    std::for_each(temporary_buffers.begin(), temporary_buffers.end(), [] (wgpu::Buffer& b) {b.destroy(); b.release();});
    temporary_bind_groups.clear();
    temporary_buffers.clear();
}

template<bool weighted_bins, bool variable_bin_width>
int WebGPUBackend<weighted_bins, variable_bin_width>::size_self_result() const {
    return self_result_count;
}

template<bool weighted_bins, bool variable_bin_width>
int WebGPUBackend<weighted_bins, variable_bin_width>::size_cross_result() const {
    return cross_result_count;
}

template class ausaxs::gpu::WebGPUBackend<false, false>;
template class ausaxs::gpu::WebGPUBackend<false, true>;
template class ausaxs::gpu::WebGPUBackend<true, false>;
template class ausaxs::gpu::WebGPUBackend<true, true>;
