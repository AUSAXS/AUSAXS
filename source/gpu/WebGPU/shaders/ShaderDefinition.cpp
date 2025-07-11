#include <gpu/WebGPU/shaders/ShaderDefinition.h>
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

ShaderDefinition::ShaderDefinition(std::string_view source, wgpu::Device device) {
    module = load_shader_module(source, device);
}