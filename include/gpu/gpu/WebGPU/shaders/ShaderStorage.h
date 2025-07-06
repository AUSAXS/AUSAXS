#pragma once

#include <gpu/WebGPU/shaders/ShaderDefinition.h>

// Storage for shaders and their layouts. The idea is to perform lazy initialization of each module, 
// since validation of the shader source code is expensive and we want to avoid it if possible. 
namespace ausaxs::gpu::shader {
    class Simple {
        public:
            Simple() = default;
            Simple(wgpu::Device device);

            ShaderDefinition weighted();
            ShaderDefinition unweighted();

            template<bool weighted_bins>
            ShaderDefinition get() {
                if constexpr (weighted_bins) {
                    return weighted();
                } else {
                    return unweighted();
                }
            }

        private:
            wgpu::Device device;
            ShaderDefinition weighted_shader;
            ShaderDefinition unweighted_shader;
    };
}