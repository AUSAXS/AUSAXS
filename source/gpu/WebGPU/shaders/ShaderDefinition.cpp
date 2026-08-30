// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/WebGPU/shaders/ShaderDefinition.h>
#include <gpu/WebGPU/GPUInstance.h>
#include <utility/Exceptions.h>
#include <utility/Logging.h>

#include <array>

using namespace ausaxs;
using namespace ausaxs::gpu;

namespace {
    wgpu::ShaderModule compile_shader_module(std::string_view source) {
        wgpu::ShaderSourceWGSL shader_source;
        shader_source.code = wgpu::StringView(source);
        shader_source.chain.next = nullptr;
        shader_source.chain.sType = wgpu::SType::ShaderSourceWGSL;

        wgpu::ShaderModuleDescriptor shader_module_descriptor{};
        shader_module_descriptor.nextInChain = &shader_source.chain;
        wgpu::ShaderModule module = GPUInstance::get().device.createShaderModule(shader_module_descriptor);
        if (!module) {throw except::runtime_error("ShaderDefinition: could not compile the shader module.");}
        logging::log("ShaderDefinition: compiled shader module.");
        return module;
    }

    // the layout is shared by all simple shaders: two atom buffers, a histogram buffer and the parameters
    wgpu::BindGroupLayout create_bind_group_layout() {
        std::array<wgpu::BindGroupLayoutEntry, 4> bindings = {wgpu::Default, wgpu::Default, wgpu::Default, wgpu::Default};

        // atom buffers
        bindings[0].binding = 0;
        bindings[0].buffer.type = wgpu::BufferBindingType::ReadOnlyStorage;
        bindings[0].visibility = wgpu::ShaderStage::Compute;

        bindings[1].binding = 1;
        bindings[1].buffer.type = wgpu::BufferBindingType::ReadOnlyStorage;
        bindings[1].visibility = wgpu::ShaderStage::Compute;

        // histogram buffer
        bindings[2].binding = 2;
        bindings[2].buffer.type = wgpu::BufferBindingType::Storage;
        bindings[2].visibility = wgpu::ShaderStage::Compute;

        // parameter buffer
        bindings[3].binding = 3;
        bindings[3].buffer.type = wgpu::BufferBindingType::Uniform;
        bindings[3].visibility = wgpu::ShaderStage::Compute;

        wgpu::BindGroupLayoutDescriptor bindgroup_layout_desc;
        bindgroup_layout_desc.entryCount = bindings.size();
        bindgroup_layout_desc.entries = bindings.data();
        return GPUInstance::get().device.createBindGroupLayout(bindgroup_layout_desc);
    }

    ShaderDefinition::ComputePipelines create_pipelines(wgpu::BindGroupLayout bind_group_layout, wgpu::ShaderModule module) {
        auto device = GPUInstance::get().device;
        auto ctr = [&] (std::string_view entry_point) {
            wgpu::PipelineLayoutDescriptor pipeline_layout_desc;
            pipeline_layout_desc.bindGroupLayoutCount = 1;
            pipeline_layout_desc.bindGroupLayouts = (WGPUBindGroupLayout*) &bind_group_layout;
            auto pipeline_layout = device.createPipelineLayout(pipeline_layout_desc);

            wgpu::ComputePipelineDescriptor pipeline_descriptor = wgpu::Default;
            pipeline_descriptor.compute.entryPoint.data = entry_point.data();
            pipeline_descriptor.compute.entryPoint.length = entry_point.size();
            pipeline_descriptor.compute.module = module;
            pipeline_descriptor.layout = pipeline_layout;
            auto pipeline = device.createComputePipeline(pipeline_descriptor);
            if (!pipeline) {throw except::runtime_error("ShaderDefinition: could not create the \"" + std::string(entry_point) + "\" pipeline.");}
            return pipeline;
        };

        return {
            .self = ctr("calculate_self"),
            .cross = ctr("calculate_cross")
        };
    }
}

ShaderDefinition::ShaderDefinition(std::string_view source) {
    module = compile_shader_module(source);
    bind_group_layout = create_bind_group_layout();
    pipelines = create_pipelines(bind_group_layout, module);
}
