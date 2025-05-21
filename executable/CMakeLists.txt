FetchContent_Declare(
	CLI11
	GIT_REPOSITORY https://github.com/CLIUtils/CLI11
	GIT_TAG 792d89286788acac125e0487f8dbde88035f7422
	GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(CLI11)

add_executable(saxs_fitter "saxs_fitter.cpp")
add_executable(em_fitter "em_fitter.cpp")
add_executable(rigidbody_optimizer "rigidbody_optimizer.cpp")

target_link_libraries(saxs_fitter PRIVATE ausaxs_core ausaxs_math CLI11::CLI11 webgpu)
target_link_libraries(em_fitter PRIVATE ausaxs_core ausaxs_math ausaxs_em CLI11::CLI11)
target_link_libraries(rigidbody_optimizer PRIVATE ausaxs_core ausaxs_math ausaxs_rigidbody CLI11::CLI11)
target_copy_webgpu_binaries(saxs_fitter)

if (GUI)
	add_subdirectory(gui)
endif()
