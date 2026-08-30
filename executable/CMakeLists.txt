# helper function to attach post-build copy of plotting scripts to an executable target.
function(add_plot_scripts_to_target _target)
	if (NOT TARGET ${_target})
		message(WARNING "add_plot_scripts_to_target: target ${_target} does not exist")
		return()
	endif()

	add_custom_command(TARGET ${_target} POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E make_directory $<TARGET_FILE_DIR:${_target}>
		COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/scripts/plot.py $<TARGET_FILE_DIR:${_target}>/plot.py
		COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/scripts/plot_helper.py $<TARGET_FILE_DIR:${_target}>/plot_helper.py
		COMMENT "Copying plotting scripts to runtime directory for ${_target}"
	)
endfunction()

add_executable(saxs_fitter "saxs_fitter.cpp")
add_executable(em_fitter "em_fitter.cpp")
add_executable(rigidbody_optimizer "rigidbody_optimizer.cpp")
add_executable(ausaxs "ausaxs.cpp")

target_link_libraries(saxs_fitter PRIVATE ausaxs_core ausaxs_api ausaxs_math ausaxs_em ausaxs_rigidbody)
target_link_libraries(em_fitter PRIVATE ausaxs_core ausaxs_api ausaxs_math ausaxs_em ausaxs_rigidbody)
target_link_libraries(rigidbody_optimizer PRIVATE ausaxs_core ausaxs_api ausaxs_math ausaxs_em ausaxs_rigidbody)
target_link_libraries(ausaxs PRIVATE ausaxs_core ausaxs_api ausaxs_math ausaxs_em ausaxs_rigidbody)

add_plot_scripts_to_target(saxs_fitter)
add_plot_scripts_to_target(em_fitter)
add_plot_scripts_to_target(rigidbody_optimizer)
add_plot_scripts_to_target(ausaxs)

if (GPU OR SYCL)
	# compares the GPU histogram kernels against the CPU one
	add_executable(gpu_debug "gpu_debug.cpp")
	target_link_libraries(gpu_debug PRIVATE ausaxs_core ausaxs_api ausaxs_math ausaxs_em ausaxs_rigidbody)
endif()

if (SYCL)
	# the SYCL kernels are a prebuilt object that has to be embedded in each final binary; see source/gpu
	foreach(_target saxs_fitter em_fitter rigidbody_optimizer ausaxs gpu_debug)
		ausaxs_add_sycl_kernels(${_target})
	endforeach()
endif()

if (GPU)
	target_copy_webgpu_binaries(gpu_debug)
	target_copy_webgpu_binaries(saxs_fitter)
	target_copy_webgpu_binaries(em_fitter)
	target_copy_webgpu_binaries(rigidbody_optimizer)
	target_copy_webgpu_binaries(ausaxs)
endif()

if (GUI)
	add_subdirectory(gui)
endif()
