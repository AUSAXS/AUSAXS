#include <gpu/WebGPU/BufferManager.h>
#include <gpu/WebGPU/BindGroups.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesData.h>
#include <constants/ConstantsAxes.h>

#include <thread>

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

wgpu::BindGroup assign(const InstanceManager& instance, wgpu::Buffer atoms, wgpu::Buffer hist) {
    assert(atoms && "Buffer is null");
    assert(hist && "Histogram buffer is null");

    std::vector<wgpu::BindGroupEntry> entries(3, wgpu::Default);
    entries[0].binding = 0;
    entries[0].buffer = atoms;
    entries[0].size = atoms.getSize();

    //? not ideal to use a "fake" buffer like this
    entries[1].binding = 1;
    entries[1].buffer = atoms;
    entries[1].size = 4;

    entries[2].binding = 2;
    entries[2].buffer = hist;
    entries[2].size = hist.getSize();

    wgpu::BindGroupDescriptor bind_group_desc;
    bind_group_desc.layout = instance.bind_group_layout;
    bind_group_desc.entryCount = entries.size();
    bind_group_desc.entries = (WGPUBindGroupEntry*) entries.data();
    return instance.device.createBindGroup(bind_group_desc);
}

wgpu::BindGroup assign(const InstanceManager& instance, wgpu::Buffer atoms1, wgpu::Buffer atoms2, wgpu::Buffer hist) {
    assert(atoms1 && "Buffer 1 is null");
    assert(atoms2 && "Buffer 2 is null");
    assert(hist && "Histogram buffer is null");

    std::vector<wgpu::BindGroupEntry> entries(3, wgpu::Default);
    entries[0].binding = 0;
    entries[0].buffer = atoms1;
    entries[0].size = atoms1.getSize();

    entries[1].binding = 1;
    entries[1].buffer = atoms2;
    entries[1].size = atoms2.getSize();

    entries[2].binding = 2;
    entries[2].buffer = hist;
    entries[2].size = hist.getSize();

    wgpu::BindGroupDescriptor bind_group_desc;
    bind_group_desc.layout = instance.bind_group_layout;
    bind_group_desc.entryCount = entries.size();
    bind_group_desc.entries = (WGPUBindGroupEntry*) entries.data();
    return instance.device.createBindGroup(bind_group_desc);
}

std::vector<wgpu::Buffer> readback_buffers;
wgpu::BindGroup BufferManager::create(wgpu::CommandEncoder encoder, const CompactCoordinates& data) {
    auto atom_buffer = create_atomic_buffer(instance.device, data);
    auto hist_buffer = create_histogram_buffer(instance.device);
    hist_buffers.emplace_back(hist_buffer);
    wgpu::BindGroup group = assign(instance, atom_buffer, hist_buffer);

    wgpu::BufferDescriptor readback_buffer_desc;
    readback_buffer_desc.size = hist_buffer.getSize();
    readback_buffer_desc.mappedAtCreation = false;
    readback_buffer_desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
    wgpu::Buffer readback_buffer = instance.device.createBuffer(readback_buffer_desc);
    encoder.copyBufferToBuffer(hist_buffer, 0, readback_buffer, 0, readback_buffer.getSize());
    readback_buffers.push_back(readback_buffer);

    return group;
}

wgpu::BindGroup BufferManager::create(wgpu::CommandEncoder encoder, const CompactCoordinates& data1, const CompactCoordinates& data2) {
    auto atom_buffer_1 = create_atomic_buffer(instance.device, data1);
    auto atom_buffer_2 = create_atomic_buffer(instance.device, data2);
    auto hist_buffer = create_histogram_buffer(instance.device);
    hist_buffers.emplace_back(hist_buffer);
    wgpu::BindGroup group = assign(instance, atom_buffer_1, atom_buffer_2, hist_buffer);

    wgpu::BufferDescriptor readback_buffer_desc;
    readback_buffer_desc.size = hist_buffer.getSize();
    readback_buffer_desc.mappedAtCreation = false;
    readback_buffer_desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
    wgpu::Buffer readback_buffer = instance.device.createBuffer(readback_buffer_desc);
    encoder.copyBufferToBuffer(hist_buffer, 0, readback_buffer, 0, readback_buffer.getSize());
    readback_buffers.push_back(readback_buffer);

    return group;
}

struct CallbackBufferData {
    wgpu::Buffer* gpu;
    std::vector<float>* cpu;
};
std::vector<double> BufferManager::merge_histograms(wgpu::Instance instance) const {
    std::vector<std::vector<float>> results;

    // //! This should be performed in parallel in a shader function
    // std::vector<wgpu::Buffer> readback_buffers;
    // {   // copy to readable location
    //     wgpu::CommandEncoder encoder = device.createCommandEncoder();
    //     for (const auto& buffer : hist_buffers) {
    //         if (!buffer) {
    //             throw std::runtime_error("BufferManager::merge_histograms: One of the histogram buffers is null.");
    //         }
    //         wgpu::BufferDescriptor readback_buffer_desc;
    //         readback_buffer_desc.size = buffer.getSize();
    //         readback_buffer_desc.mappedAtCreation = false;
    //         readback_buffer_desc.usage = wgpu::BufferUsage::CopyDst | wgpu::BufferUsage::MapRead;
    //         wgpu::Buffer readback_buffer = device.createBuffer(readback_buffer_desc);
            
    //         encoder.copyBufferToBuffer(buffer, 0, readback_buffer, 0, readback_buffer.getSize());
    //         readback_buffers.emplace_back(std::move(readback_buffer));
    //     }

    //     wgpu::CommandBufferDescriptor command_buffer_desc;
    //     command_buffer_desc.nextInChain = nullptr;
    //     auto command_buffer = encoder.finish(command_buffer_desc);
    //     encoder.release();

    //     device.getQueue().submit(command_buffer);
    //     command_buffer.release();
    // }

    {   // copy to cpu
        for (auto& buffer : readback_buffers) {
            std::vector<float> result(constants::axes::d_vals.size());
            CallbackBufferData callback_data{&buffer, &result};

            bool done = false;
            wgpu::BufferMapCallbackInfo map_callback;
            map_callback.mode = wgpu::CallbackMode::AllowProcessEvents;
            map_callback.userdata1 = &callback_data;
            map_callback.userdata2 = &done;
            map_callback.callback = [] (WGPUMapAsyncStatus, WGPUStringView, void* p_readback_buffer, void* p_done) -> void {
                CallbackBufferData buffers = *reinterpret_cast<CallbackBufferData*>(p_readback_buffer);
                bool* done = reinterpret_cast<bool*>(p_done);
                const float* output = (const float*) buffers.gpu->getConstMappedRange(0, buffers.gpu->getSize());
                buffers.cpu->assign(output, output + buffers.cpu->size());
                buffers.gpu->unmap();
                *done = true;
            };
            buffer.mapAsync(wgpu::MapMode::Read, 0, buffer.getSize(), map_callback);
            instance.processEvents();

            int N = 0;
            while(!done && ++N < 100) {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                instance.processEvents();                
            }
            if (!done) {
                throw std::runtime_error("BufferManager::merge_histograms: Readback buffer did not finish in time.");
            }

            results.emplace_back(std::move(result));
            buffer.release();
        }
    }

    // finally, merge the results
    std::vector<double> result(constants::axes::d_vals.size());
    for (const auto& res : results) {
        std::transform(result.begin(), result.end(), res.begin(), result.begin(), std::plus<double>());
    }
    return result;
}