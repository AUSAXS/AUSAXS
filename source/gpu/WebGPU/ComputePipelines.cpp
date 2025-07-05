#include "gpu/WebGPU/InstanceManager.h"
#include <gpu/WebGPU/ComputePipelines.h>
#include <gpu/WebGPU/BindGroups.h>
#include <io/ExistingFile.h>
#include <utility/Logging.h>

#include <fstream>
#include <iterator>
#include <stdexcept>

using namespace ausaxs;
using namespace ausaxs::gpu;

wgpu::ShaderModule load_shader_module(const io::ExistingFile& path, wgpu::Device device) {
    // read whole source file into a string
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("load_shader_module: Could not open shader file: \"" + path.path() + "\"");
    }
    std::string source(std::istreambuf_iterator<char>{file}, {});

    // define shader type
    wgpu::ShaderSourceWGSL shader_source;
    shader_source.code = wgpu::StringView(source);
    shader_source.chain.next = nullptr;
    shader_source.chain.sType = wgpu::SType::ShaderSourceWGSL;

    // create shader module descriptor
    wgpu::ShaderModuleDescriptor shader_module_descriptor{};
    shader_module_descriptor.nextInChain = &shader_source.chain;
    wgpu::ShaderModule module = device.createShaderModule(shader_module_descriptor);
    logging::log("WebGPU: Succesfully loaded shader module.");
    return module;
}

ComputePipelines::Pipelines ComputePipelines::create(const InstanceManager& instance) {
    static wgpu::ShaderModule compute_shader_module = load_shader_module("executable/gpu_debug.wgsl", instance.device);

    auto ctr = [&] (std::string_view shader_name) {
        wgpu::PipelineLayoutDescriptor pipeline_layout_desc;
        pipeline_layout_desc.bindGroupLayoutCount = 1;
        pipeline_layout_desc.bindGroupLayouts = (WGPUBindGroupLayout*) &instance.bind_group_layout;
        auto pipeline_layout = instance.device.createPipelineLayout(pipeline_layout_desc);

        wgpu::ComputePipelineDescriptor pipeline_descriptor = wgpu::Default;
        pipeline_descriptor.compute.entryPoint.data = shader_name.data();
        pipeline_descriptor.compute.entryPoint.length = shader_name.size();
        pipeline_descriptor.compute.module = compute_shader_module;
        pipeline_descriptor.layout = pipeline_layout;
        return instance.device.createComputePipeline(pipeline_descriptor);
    };

    return {
        .self=ctr("calculate_self"),
        .cross=ctr("calculate_cross")
    };
}