#include <gpu/WebGPU/shaders/ShaderStorage.h>

using namespace ausaxs::gpu;

shader::Simple::Simple(wgpu::Device device) : device(device) {};

wgpu::BindGroupLayoutDescriptor simple_layout() {
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

    static wgpu::BindGroupLayoutDescriptor bindgroup_layout_desc = [] () {
        wgpu::BindGroupLayoutDescriptor bindgroup_layout_desc;
        bindgroup_layout_desc.entryCount = bindings.size();
        bindgroup_layout_desc.entries = bindings.data();
        return bindgroup_layout_desc;
    }();

    return bindgroup_layout_desc;
}

ShaderDefinition shader::Simple::weighted() {
    unweighted_shader = ShaderDefinition("source/gpu/WebGPU/shaders/simple/WebGPUSimpleWeighted.wgsl", device);
    unweighted_shader.layout = device.createBindGroupLayout(simple_layout());
    return unweighted_shader;
}

ShaderDefinition shader::Simple::unweighted() {
    weighted_shader = ShaderDefinition("source/gpu/WebGPU/shaders/simple/WebGPUSimpleUnweighted.wgsl", device);
    weighted_shader.layout = device.createBindGroupLayout(simple_layout());
    return weighted_shader;
}