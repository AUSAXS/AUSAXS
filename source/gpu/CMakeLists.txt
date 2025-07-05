add_library(ausaxs_gpu OBJECT 
	"WebGPUSimple.cpp"

	"WebGPU/BindGroups.cpp"
	"WebGPU/Buffers.cpp"
	"WebGPU/ComputePipelines.cpp"
	"WebGPU/InstanceManager.cpp"
)

target_include_directories(ausaxs_gpu PUBLIC "${CMAKE_SOURCE_DIR}/include/gpu")
target_link_libraries(ausaxs_gpu ausaxs_core webgpu)
set_target_properties(ausaxs_gpu PROPERTIES POSITION_INDEPENDENT_CODE 1)
