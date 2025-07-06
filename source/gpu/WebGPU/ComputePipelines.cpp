#include <gpu/WebGPU/ComputePipelines.h>
#include <gpu/WebGPU/InstanceManager.h>
#include <gpu/WebGPU/BindGroups.h>
#include <io/ExistingFile.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::gpu;

template<bool weighted_bins>
ComputePipelines<weighted_bins>::Pipelines ComputePipelines<weighted_bins>::create(wgpu::Device device, ShaderDefinition shader) {
    auto ctr = [&] (std::string_view shader_name) {
        wgpu::PipelineLayoutDescriptor pipeline_layout_desc;
        pipeline_layout_desc.bindGroupLayoutCount = 1;
        pipeline_layout_desc.bindGroupLayouts = (WGPUBindGroupLayout*) &shader.layout;
        auto pipeline_layout = device.createPipelineLayout(pipeline_layout_desc);

        wgpu::ComputePipelineDescriptor pipeline_descriptor = wgpu::Default;
        pipeline_descriptor.compute.entryPoint.data = shader_name.data();
        pipeline_descriptor.compute.entryPoint.length = shader_name.size();
        pipeline_descriptor.compute.module = shader.module;
        pipeline_descriptor.layout = pipeline_layout;
        return device.createComputePipeline(pipeline_descriptor);
    };

    return {
        .self = ctr("calculate_self"),
        .cross = ctr("calculate_cross")
    };
}

template struct ausaxs::gpu::ComputePipelines<false>;
template struct ausaxs::gpu::ComputePipelines<true>;