// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/WebGPU/Buffers.h>
#include <gpu/WebGPU/GPUInstance.h>
#include <settings/HistogramSettings.h>

using namespace ausaxs;
using namespace ausaxs::gpu;
using namespace ausaxs::hist::detail;

template<bool weighted_bins> template<bool variable_bin_width>
wgpu::Buffer Buffers<weighted_bins>::create_atom_buffer(const CompactCoordinates<variable_bin_width>& atoms) {
    using AtomType = std::remove_cvref_t<decltype(atoms[0])>;
    static_assert(sizeof(AtomType) == 4*sizeof(float), "The atom type must be 4 floats in size to match the shader layout.");

    auto& gpu_instance = GPUInstance::get();
    wgpu::BufferDescriptor atom_buffer_desc;
    atom_buffer_desc.size = atoms.size()*sizeof(AtomType);
    atom_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
    wgpu::Buffer atom_buffer = gpu_instance.device.createBuffer(atom_buffer_desc);
    gpu_instance.queue.writeBuffer(atom_buffer, 0, atoms.get_data().data(), atom_buffer.getSize());
    return atom_buffer;
}

template<bool weighted_bins>
wgpu::Buffer Buffers<weighted_bins>::dummy_buffer() {
    static wgpu::Buffer dummy = [] () {
        wgpu::BufferDescriptor dummy_buffer_desc;
        dummy_buffer_desc.size = 4*sizeof(float); // a single, unused atom
        dummy_buffer_desc.usage = wgpu::BufferUsage::Storage;
        return GPUInstance::get().device.createBuffer(dummy_buffer_desc);
    }();
    return dummy;
}

template<bool weighted_bins>
wgpu::Buffer Buffers<weighted_bins>::create_histogram_buffer() {
    wgpu::BufferDescriptor histogram_buffer_desc;
    histogram_buffer_desc.size = settings::axes::bin_count*sizeof(HistogramType);
    histogram_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
    return GPUInstance::get().device.createBuffer(histogram_buffer_desc);
}

template<bool weighted_bins>
wgpu::Buffer Buffers<weighted_bins>::create_readback_buffer() {
    wgpu::BufferDescriptor readback_buffer_desc;
    readback_buffer_desc.size = settings::axes::bin_count*sizeof(HistogramType);
    readback_buffer_desc.mappedAtCreation = false;
    readback_buffer_desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
    return GPUInstance::get().device.createBuffer(readback_buffer_desc);
}

template<bool weighted_bins>
wgpu::Buffer Buffers<weighted_bins>::create_parameter_buffer(const ShaderParameters& parameters) {
    auto& gpu_instance = GPUInstance::get();
    wgpu::BufferDescriptor parameter_buffer_desc;
    parameter_buffer_desc.size = sizeof(ShaderParameters);
    parameter_buffer_desc.usage = wgpu::BufferUsage::Uniform | wgpu::BufferUsage::CopyDst;
    wgpu::Buffer parameter_buffer = gpu_instance.device.createBuffer(parameter_buffer_desc);
    gpu_instance.queue.writeBuffer(parameter_buffer, 0, &parameters, sizeof(ShaderParameters));
    return parameter_buffer;
}

template struct ausaxs::gpu::Buffers<false>;
template struct ausaxs::gpu::Buffers<true>;
template wgpu::Buffer ausaxs::gpu::Buffers<false>::create_atom_buffer<false>(const CompactCoordinates<false>&);
template wgpu::Buffer ausaxs::gpu::Buffers<false>::create_atom_buffer<true>(const CompactCoordinates<true>&);
template wgpu::Buffer ausaxs::gpu::Buffers<true>::create_atom_buffer<false>(const CompactCoordinates<false>&);
template wgpu::Buffer ausaxs::gpu::Buffers<true>::create_atom_buffer<true>(const CompactCoordinates<true>&);
