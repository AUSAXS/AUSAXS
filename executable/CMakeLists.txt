FetchContent_Declare(
	CLI11
	GIT_REPOSITORY https://github.com/CLIUtils/CLI11
	GIT_TAG 792d89286788acac125e0487f8dbde88035f7422
	GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(CLI11)

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

target_link_libraries(saxs_fitter PRIVATE ausaxs_core ausaxs_math CLI11::CLI11)
target_link_libraries(em_fitter PRIVATE ausaxs_core ausaxs_math ausaxs_em CLI11::CLI11)
target_link_libraries(rigidbody_optimizer PRIVATE ausaxs_core ausaxs_math ausaxs_rigidbody CLI11::CLI11)

add_plot_scripts_to_target(saxs_fitter)
add_plot_scripts_to_target(em_fitter)
add_plot_scripts_to_target(rigidbody_optimizer)

if (GUI)
	add_subdirectory(gui)
endif()