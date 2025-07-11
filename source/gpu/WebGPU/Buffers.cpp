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
    atom_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst;
    wgpu::Buffer atom_buffer = device.createBuffer(atom_buffer_desc);
    queue.writeBuffer(atom_buffer, 0, atoms.get_data().data(), atom_buffer.getSize());
    return atom_buffer;
}

wgpu::Buffer create_dummy_buffer(wgpu::Device device) {
    wgpu::BufferDescriptor dummy_buffer_desc;
    dummy_buffer_desc.size = 4; // smallest possible buffer size
    dummy_buffer_desc.usage = wgpu::BufferUsage::Storage;
    return device.createBuffer(dummy_buffer_desc);
}

template<bool weighted_bins>
wgpu::Buffer create_histogram_buffer(wgpu::Device device) {
    wgpu::BufferDescriptor histogram_buffer_desc;
    histogram_buffer_desc.size = constants::axes::d_axis.bins*sizeof(typename Buffers<weighted_bins>::HistogramType);
    histogram_buffer_desc.usage = wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopySrc;
    return device.createBuffer(histogram_buffer_desc);
}

template<bool weighted_bins>
BufferInstance Buffers<weighted_bins>::create(wgpu::Device device, const CompactCoordinates& data) {
    auto atom_buffer = create_atomic_buffer(device, data);
    auto dummy_buffer = create_dummy_buffer(device);
    auto hist_buffer = create_histogram_buffer<weighted_bins>(device);
    return {atom_buffer, dummy_buffer, hist_buffer};
}

template<bool weighted_bins>
BufferInstance Buffers<weighted_bins>::create(wgpu::Device device, const CompactCoordinates& data1, const CompactCoordinates& data2) {
    auto atomic_1 = create_atomic_buffer(device, data1);
    auto atomic_2 = create_atomic_buffer(device, data2);
    auto hist = create_histogram_buffer<weighted_bins>(device);
    return {atomic_1, atomic_2, hist};
}

template struct ausaxs::gpu::Buffers<false>;
template struct ausaxs::gpu::Buffers<true>;