# The GPU sources are compiled into ausaxs_core rather than a separate library, since the histogram
# managers pick their kernel at runtime and would otherwise create a circular dependency between the two.

target_include_directories(ausaxs_core PUBLIC "${CMAKE_SOURCE_DIR}/include/gpu")

if (GPU)
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

	target_link_libraries(ausaxs_core webgpu)
	target_compile_definitions(ausaxs_core PUBLIC AUSAXS_GPU)
endif()

# The SYCL kernels have to be built by the SYCL compiler rather than the one used for the rest of the
# project, so they are compiled on their own and the resulting object is linked into ausaxs_core. This
# is why everything they exchange with the rest of the library is trivially copyable; see SYCLBackend.h.
if (SYCL)
	find_program(ACPP_EXECUTABLE
		NAMES acpp syclcc
		HINTS ${ACPP_ROOT}/bin $ENV{ACPP_ROOT}/bin /opt/AdaptiveCpp/bin /usr/local/bin
		DOC "The AdaptiveCpp compiler driver"
	)
	if (NOT ACPP_EXECUTABLE)
		message(FATAL_ERROR "SYCL support was requested, but no acpp was found. Set -DACPP_ROOT=<AdaptiveCpp install>.")
	endif()
	get_filename_component(_acpp_bin "${ACPP_EXECUTABLE}" DIRECTORY)
	get_filename_component(_acpp_prefix "${_acpp_bin}" DIRECTORY)
	set(_acpp_lib "${_acpp_prefix}/lib")

	# acpp records the clang it was configured with, but its plugin can only be loaded by the clang
	# matching the LLVM it was built against. These drift apart when the system compiler is upgraded,
	# so fall back to the versioned clang of that LLVM if the recorded one no longer fits.
	set(_acpp_clang "${ACPP_CLANG}")
	if (NOT _acpp_clang AND EXISTS "${_acpp_prefix}/etc/AdaptiveCpp/acpp-core.json")
		file(READ "${_acpp_prefix}/etc/AdaptiveCpp/acpp-core.json" _acpp_config)
		string(JSON _plugin_llvm ERROR_VARIABLE _ignored GET "${_acpp_config}" "plugin-llvm-version-major")
		string(JSON _configured_clang ERROR_VARIABLE _ignored GET "${_acpp_config}" "default-clang")
		if (_plugin_llvm AND _configured_clang)
			execute_process(
				COMMAND "${_configured_clang}" -dumpversion
				OUTPUT_VARIABLE _clang_version OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
			)
			if (NOT _clang_version MATCHES "^${_plugin_llvm}\\.")
				find_program(_versioned_clang "clang++-${_plugin_llvm}")
				if (_versioned_clang)
					set(_acpp_clang "${_versioned_clang}")
					message(STATUS
						"AdaptiveCpp was configured with ${_configured_clang} (${_clang_version}), which no longer "
						"matches its LLVM ${_plugin_llvm} plugin; using ${_acpp_clang} instead."
					)
				else()
					message(WARNING
						"AdaptiveCpp needs clang ${_plugin_llvm} to load its plugin, but ${_configured_clang} is "
						"${_clang_version} and no clang++-${_plugin_llvm} was found. Set -DACPP_CLANG=<clang>."
					)
				endif()
			endif()
		endif()
	endif()
	if (_acpp_clang)
		list(APPEND _acpp_flags "--acpp-clang=${_acpp_clang}")
	endif()

	set(_sycl_source "${CMAKE_CURRENT_SOURCE_DIR}/SYCL/SYCLBackend.cpp")
	set(_sycl_object "${CMAKE_CURRENT_BINARY_DIR}/SYCLBackend.o")
	set(_acpp_command ${ACPP_EXECUTABLE}
		--acpp-targets=${ACPP_TARGETS} ${_acpp_flags}
		-c -O3 -fPIC -std=c++20
		-I "${CMAKE_SOURCE_DIR}/include/gpu"
		-o ${_sycl_object} ${_sycl_source}
	)

	# compile once at configure time so the object exists when ausaxs_core is defined; the custom command
	# below keeps it up to date afterwards. this is needed because the target lives in another directory,
	# where the GENERATED source property is not visible under the policies of our cmake_minimum_required.
	if (NOT EXISTS ${_sycl_object})
		message(STATUS "Compiling the SYCL kernels for ${ACPP_TARGETS}")
		execute_process(COMMAND ${_acpp_command} COMMAND_ERROR_IS_FATAL ANY)
	endif()

	add_custom_command(
		OUTPUT ${_sycl_object}
		COMMAND ${_acpp_command}
		DEPENDS ${_sycl_source} "${CMAKE_SOURCE_DIR}/include/gpu/gpu/SYCL/SYCLBackend.h"
		COMMENT "Compiling the SYCL kernels for ${ACPP_TARGETS}"
		VERBATIM
	)
	add_custom_target(ausaxs_sycl_kernels DEPENDS ${_sycl_object})
	add_dependencies(ausaxs_core ausaxs_sycl_kernels)

	# ausaxs_core is an object library, and objects that it did not compile itself are not carried by
	# $<TARGET_OBJECTS:ausaxs_core>. Passing the object on as a link item of ausaxs_core would carry it,
	# but it would also propagate to everything downstream: a target linking the shared libausaxs, which
	# already embeds it, would receive a second copy of the device code. AdaptiveCpp registers both under
	# the same id and reports "hcf object id collision [...] some kernels might be unavailable". The
	# object is therefore embedded in each final binary as a private source instead, exactly once each.
	set(AUSAXS_SYCL_OBJECT "${_sycl_object}" CACHE INTERNAL "The compiled SYCL kernel object")

	# for the binaries defined after this directory; those defined before it are handled just below.
	# a target that links ausaxs_core without calling this fails to link, so it cannot pass unnoticed.
	function(ausaxs_add_sycl_kernels _target)
		target_sources(${_target} PRIVATE "${AUSAXS_SYCL_OBJECT}")
		add_dependencies(${_target} ausaxs_sycl_kernels)
	endfunction()

	ausaxs_add_sycl_kernels(libausaxs)
	ausaxs_add_sycl_kernels(libausaxs_static)

	# the acpp runtime is shared, so there is no harm in this propagating the usual way
	target_link_libraries(ausaxs_core "-L${_acpp_lib}" acpp-rt acpp-common "-Wl,-rpath,${_acpp_lib}")
	target_compile_definitions(ausaxs_core PUBLIC AUSAXS_GPU_SYCL)
endif()
