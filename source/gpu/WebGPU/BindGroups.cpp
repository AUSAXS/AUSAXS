#include <gpu/WebGPU/BindGroups.h>
#include <utility/Logging.h>

using namespace ausaxs::gpu;

wgpu::BindGroupLayout BindGroups::get(wgpu::Device device) {
    static wgpu::Device ptr = device;
    if (device != ptr) {throw std::runtime_error("ComputePipelines::create: Device pointer has been changed.");}

    // layout is fixed for the current shader, so we can use a static array
    static std::array<wgpu::BindGroupLayoutEntry, 3> bindings = [] () {
        std::array<wgpu::BindGroupLayoutEntry, 3> b = {wgpu::Default, wgpu::Default, wgpu::Default};

        // atom buffers
        b[0].binding = 0;
        b[0].buffer.type = wgpu::BufferBindingType::ReadOnlyStorage;
        b[0].visibility = wgpu::ShaderStage::Compute;

        b[1].binding = 1;
        b[1].buffer.type = wgpu::BufferBindingType::ReadOnlyStorage;
        b[1].visibility = wgpu::ShaderStage::Compute;

        // histogram buffer
        b[2].binding = 2;
        b[2].buffer.type = wgpu::BufferBindingType::Storage;
        b[2].visibility = wgpu::ShaderStage::Compute;
        return b;
    }();

    // bindgroup layout can also be static
    static wgpu::BindGroupLayoutDescriptor bindgroup_layout = [] () {
        wgpu::BindGroupLayoutDescriptor bindgroup_layout_desc;
        bindgroup_layout_desc.entryCount = bindings.size();
        bindgroup_layout_desc.entries = bindings.data();
        return bindgroup_layout_desc;
    }();

    // create the bind group
    static wgpu::BindGroupLayout bind_group = ptr.createBindGroupLayout(bindgroup_layout);
    logging::log("WebGPU: Created bind group layout.");
    return bind_group;
}