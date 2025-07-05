#include <gpu/WebGPU/Buffers.h>

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesData.h>

using namespace ausaxs;
using namespace ausaxs::gpu;
using namespace ausaxs::hist::detail;

wgpu::Buffer create_atomic_buffer(wgpu::Device device, const CompactCoordinates& atoms) {
    static_assert(sizeof(CompactCoordinatesData) == 4*sizeof(float), "CompactCoordinatesData must be 4 floats in size.");

    auto queue = device.getQueue();
    wgpu::BufferDescriptor atom_buffer_desc;
    atom_buffer_desc.size = atoms.size()*sizeof(CompactCoordinatesData);
    atom_buffer_desc.mappedAtCreation = false;
    atom_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
    wgpu::Buffer atom_buffer = device.createBuffer(atom_buffer_desc);
    queue.writeBuffer(atom_buffer, 0, atoms.get_data().data(), atom_buffer.getSize());
    return atom_buffer;
}

wgpu::Buffer create_histogram_buffer(wgpu::Device device) {
    wgpu::BufferDescriptor histogram_buffer_desc;
    histogram_buffer_desc.size = constants::axes::d_vals.size()*sizeof(float);
    histogram_buffer_desc.mappedAtCreation = false;
    histogram_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
    return device.createBuffer(histogram_buffer_desc);
}

Buffers::BufferInstance Buffers::create(wgpu::Device device, const CompactCoordinates& data) {
    auto atom_buffer = create_atomic_buffer(device, data);
    auto hist_buffer = create_histogram_buffer(device);
    return {atom_buffer, {}, hist_buffer};
}

Buffers::BufferInstance Buffers::create(wgpu::Device device, const CompactCoordinates& data1, const CompactCoordinates& data2) {
    auto atomic_1 = create_atomic_buffer(device, data1);
    auto atomic_2 = create_atomic_buffer(device, data2);
    auto hist = create_histogram_buffer(device);
    return {atomic_1, atomic_2, hist};
}
