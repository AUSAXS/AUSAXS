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

ShaderDefinition::ComputePipelines simple_pipeline(wgpu::Device device, wgpu::BindGroupLayout bind_group_layout, wgpu::ShaderModule module) {
    auto ctr = [&] (std::string_view shader_name) {
        wgpu::PipelineLayoutDescriptor pipeline_layout_desc;
        pipeline_layout_desc.bindGroupLayoutCount = 1;
        pipeline_layout_desc.bindGroupLayouts = (WGPUBindGroupLayout*) &bind_group_layout;
        auto pipeline_layout = device.createPipelineLayout(pipeline_layout_desc);

        wgpu::ComputePipelineDescriptor pipeline_descriptor = wgpu::Default;
        pipeline_descriptor.compute.entryPoint.data = shader_name.data();
        pipeline_descriptor.compute.entryPoint.length = shader_name.size();
        pipeline_descriptor.compute.module = module;
        pipeline_descriptor.layout = pipeline_layout;
        return device.createComputePipeline(pipeline_descriptor);
    };

    return {
        .self = ctr("calculate_self"),
        .cross = ctr("calculate_cross")
    };
}

ShaderDefinition shader::Simple::weighted() {
    unweighted_shader = ShaderDefinition("source/gpu/WebGPU/shaders/simple/WebGPUSimpleWeighted.wgsl", device);
    unweighted_shader.bind_group_layout = device.createBindGroupLayout(simple_layout());
    unweighted_shader.pipelines = simple_pipeline(device, unweighted_shader.bind_group_layout, unweighted_shader.module);
    return unweighted_shader;
}

ShaderDefinition shader::Simple::unweighted() {
    weighted_shader = ShaderDefinition("source/gpu/WebGPU/shaders/simple/WebGPUSimpleUnweighted.wgsl", device);
    weighted_shader.bind_group_layout = device.createBindGroupLayout(simple_layout());
    weighted_shader.pipelines = simple_pipeline(device, weighted_shader.bind_group_layout, weighted_shader.module);
    return weighted_shader;
}