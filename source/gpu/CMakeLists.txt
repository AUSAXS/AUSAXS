# The GPU sources are compiled into ausaxs_core rather than a separate library, since the histogram
# managers pick their kernel at runtime and would otherwise create a circular dependency between the two.

set(_shader_dir "${CMAKE_CURRENT_SOURCE_DIR}/WebGPU/shaders/simple")
set(_shader_names "simple_unweighted;simple_weighted")
set(_shader_files "${_shader_dir}/WebGPUSimpleUnweighted.wgsl;${_shader_dir}/WebGPUSimpleWeighted.wgsl")
set(_generated_shaders "${CMAKE_CURRENT_BINARY_DIR}/generated/EmbeddedShaders.cpp")

# generate once at configure time so the file exists when ausaxs_core is defined; the custom command
# below keeps it up to date afterwards. this is needed because the target lives in another directory,
# where the GENERATED source property is not visible under the policies of our cmake_minimum_required.
execute_process(COMMAND ${CMAKE_COMMAND}
	-DOUTPUT=${_generated_shaders}
	"-DNAMES=${_shader_names}"
	"-DFILES=${_shader_files}"
	-P "${CMAKE_CURRENT_SOURCE_DIR}/embed_shaders.cmake"
	COMMAND_ERROR_IS_FATAL ANY
)

add_custom_command(
	OUTPUT ${_generated_shaders}
	COMMAND ${CMAKE_COMMAND}
		-DOUTPUT=${_generated_shaders}
		"-DNAMES=${_shader_names}"
		"-DFILES=${_shader_files}"
		-P "${CMAKE_CURRENT_SOURCE_DIR}/embed_shaders.cmake"
	DEPENDS ${_shader_files} "${CMAKE_CURRENT_SOURCE_DIR}/embed_shaders.cmake"
	COMMENT "Embedding wgsl shaders"
	VERBATIM
)
# the generated file is consumed by a target defined in another directory, so it must be marked explicitly
set_source_files_properties(${_generated_shaders} PROPERTIES GENERATED TRUE)
add_custom_target(ausaxs_shaders DEPENDS ${_generated_shaders})
add_dependencies(ausaxs_core ausaxs_shaders)

target_sources(ausaxs_core PRIVATE
	${_generated_shaders}

	"WebGPU/BufferManager.cpp"
	"WebGPU/Buffers.cpp"
	"WebGPU/GPUInstance.cpp"
	"WebGPU/WebGPU.cpp"
	"WebGPU/WebGPUBackend.cpp"

	"WebGPU/shaders/ShaderDefinition.cpp"
	"WebGPU/shaders/ShaderStorage.cpp"
)

target_include_directories(ausaxs_core PUBLIC "${CMAKE_SOURCE_DIR}/include/gpu")
target_link_libraries(ausaxs_core webgpu)
target_compile_definitions(ausaxs_core PUBLIC AUSAXS_GPU)
