#include <gpu/WebGPU/WebGPU.h>
#include <gpu/WebGPU/GPUInstance.h>
#include <gpu/WebGPU/Buffers.h>
#include <gpu/WebGPU/BufferManager.h>
#include <gpu/WebGPU/shaders/ShaderStorage.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/SimpleKernel.h>

using namespace ausaxs;
using namespace ausaxs::data;

namespace ausaxs::gpu {
    template<bool weighted_bins>
    class WebGPUBackend {
        public:
            WebGPUBackend();

            int submit_self(const hist::detail::CompactCoordinates& atoms, int merge_id = -1);
            int submit_cross(const hist::detail::CompactCoordinates& atom_1, const hist::detail::CompactCoordinates& atom_2, int merge_id = -1);
            hist::distance_calculator::SimpleKernel<weighted_bins>::run_result run();

        private:
            GPUInstance instance;
            BufferManager<weighted_bins> buffer_manager;
            std::vector<wgpu::Buffer> garbage_collector; // track buffer ptrs to be released after run
            shader::Simple shaders;

            void initialize();
            wgpu::BindGroup assign_buffers(wgpu::Device device, wgpu::Buffer buffer1, wgpu::Buffer buffer2, wgpu::Buffer histogram_buffer);
    };
}