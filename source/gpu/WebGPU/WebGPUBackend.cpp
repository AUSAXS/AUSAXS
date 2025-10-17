#include <gpu/WebGPU/WebGPUBackend.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::gpu;

template<bool weighted_bins>
WebGPUBackend<weighted_bins>::WebGPUBackend() {initialize();}

template<bool weighted_bins>
void WebGPUBackend<weighted_bins>::initialize() {
    shaders = shader::Simple(instance.device);
}

template<bool weighted_bins>
hist::distance_calculator::SimpleKernel<weighted_bins>::run_result WebGPUBackend<weighted_bins>::run() {
    auto res = buffer_manager.merge(instance);
    std::for_each(garbage_collector.begin(), garbage_collector.end(), [](wgpu::Buffer& buffer) {buffer.release();});
    garbage_collector.clear();
    return buffer_manager.merge(instance);
}

template<bool weighted_bins>
wgpu::BindGroup WebGPUBackend<weighted_bins>::assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer) {
    assert(buffer1 && "Buffer 1 is null");
    assert(buffer2 && "Buffer 2 is null");
    assert(histogram_buffer && "Histogram buffer is null");
    std::cout << "Assigning buffers to bind group..." << std::endl;

    std::vector<wgpu::BindGroupEntry> entries(3, wgpu::Default);
    entries[0].binding = 0;
    entries[0].buffer = buffer1;
    entries[0].size = buffer1.getSize();

    entries[1].binding = 1;
    entries[1].buffer = buffer2;
    entries[1].size = buffer2.getSize();

    entries[2].binding = 2;
    entries[2].buffer = histogram_buffer;
    entries[2].size = histogram_buffer.getSize();

    wgpu::BindGroupDescriptor bind_group_desc;
    bind_group_desc.layout = shaders.get<weighted_bins>().bind_group_layout;
    bind_group_desc.entryCount = entries.size();
    bind_group_desc.entries = entries.data();
    return device.createBindGroup(bind_group_desc);
}

void submit_to_gpu(wgpu::Device device, wgpu::BindGroup bind_group, wgpu::ComputePipeline pipeline, unsigned int atom_count) {
    int workgroups = std::ceil(static_cast<double>(atom_count)/64);

    wgpu::CommandEncoder encoder = device.createCommandEncoder();
    auto compute_pass = encoder.beginComputePass(wgpu::Default);
    compute_pass.setPipeline(pipeline);
    compute_pass.setBindGroup(0, bind_group, 0, nullptr);
    compute_pass.dispatchWorkgroups(workgroups, 1, 1);
    compute_pass.end();

    auto command_buffer = encoder.finish();
    device.getQueue().submit(command_buffer);
    command_buffer.release();
    encoder.release();
    logging::log("WebGPU: Submitted " + std::to_string(workgroups) + " workgroups of 64 threads each to the GPU.");
}

void add_to_garbage_collector(std::vector<wgpu::Buffer>& garbage_collector, const BufferInstance& buffers) {
    garbage_collector.emplace_back(buffers.atomic_1);
    garbage_collector.emplace_back(buffers.atomic_2);
    garbage_collector.emplace_back(buffers.histogram);
}

template<bool weighted_bins>
int WebGPUBackend<weighted_bins>::submit_self(const hist::detail::CompactCoordinates& atoms, int merge_id) {
    auto buffers = Buffers<weighted_bins>::create(instance.device, atoms);
    int index = buffer_manager.manage_self(buffers.histogram, merge_id);
    add_to_garbage_collector(garbage_collector, buffers);

    wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atomic_1, buffers.atomic_2, buffers.histogram);
    submit_to_gpu(instance.device, bind_group, shaders.get<weighted_bins>().pipelines.self, atoms.size());
    return index;
}

template<bool weighted_bins>
int WebGPUBackend<weighted_bins>::submit_cross(const hist::detail::CompactCoordinates& atoms1, const hist::detail::CompactCoordinates& atoms2, int merge_id) {
    auto buffers = Buffers<weighted_bins>::create(instance.device, atoms1, atoms2);
    int index = buffer_manager.manage_cross(buffers.histogram, merge_id);
    add_to_garbage_collector(garbage_collector, buffers);

    wgpu::BindGroup bind_group = assign_buffers(instance.device, buffers.atomic_1, buffers.atomic_2, buffers.histogram);
    submit_to_gpu(instance.device, bind_group, shaders.get<weighted_bins>().pipelines.cross, atoms1.size());
    return index;
}

template class ausaxs::gpu::WebGPUBackend<false>;
template class ausaxs::gpu::WebGPUBackend<true>;