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
add_executable(bin_size_analysis "bin_size_analysis.cpp")
add_executable(exv_comparison "exv_comparison.cpp")
add_executable(fit_all_exv "fit_all_exv.cpp")
add_executable(inspect_map "inspect_map.cpp")
add_executable(plot_profiles "plot_profiles.cpp")
add_executable(thea "thea.cpp")
add_executable(vary_exv_radius "vary_exv_radius.cpp")
add_executable(vary_grid_radii "vary_grid_radii.cpp")
add_executable(pdb_to_single_atom_type "pdb_to_single_atom_type.cpp")

target_link_libraries(saxs_fitter PRIVATE ausaxs_core ausaxs_math CLI11::CLI11)
target_link_libraries(em_fitter PRIVATE ausaxs_core ausaxs_math ausaxs_em CLI11::CLI11)
target_link_libraries(rigidbody_optimizer PRIVATE ausaxs_core ausaxs_math ausaxs_rigidbody CLI11::CLI11)
target_link_libraries(exv_comparison PRIVATE ausaxs_core ausaxs_math CLI11::CLI11)
target_link_libraries(fit_all_exv PRIVATE ausaxs_core ausaxs_math CLI11::CLI11)
target_link_libraries(inspect_map PRIVATE ausaxs_core ausaxs_math CLI11::CLI11)
target_link_libraries(plot_profiles PRIVATE ausaxs_core ausaxs_math CLI11::CLI11)
target_link_libraries(thea PRIVATE ausaxs_core ausaxs_math CLI11::CLI11)
target_link_libraries(vary_exv_radius PRIVATE ausaxs_core ausaxs_math CLI11::CLI11)
target_link_libraries(vary_grid_radii PRIVATE ausaxs_core ausaxs_math CLI11::CLI11)
target_link_libraries(pdb_to_single_atom_type PRIVATE ausaxs_core ausaxs_math)

add_subdirectory(md)
if (GUI)
	add_subdirectory(gui)
endif()
