add_library(ausaxs_api OBJECT
    "ObjectStorage.cpp"

    "api_helper.cpp"
    "api_sasview.cpp"

    "pyausaxs/api_data.cpp"
    "pyausaxs/api_fit.cpp"
    "pyausaxs/api_form_factor.cpp"
    "pyausaxs/api_io.cpp"
    "pyausaxs/api_iterative_fit.cpp"
    "pyausaxs/api_molecule.cpp"
    "pyausaxs/api_output.cpp"
    "pyausaxs/api_pdb.cpp"
    "pyausaxs/api_rigidbody.cpp"
    "pyausaxs/api_settings.cpp"
    "pyausaxs/api_shape.cpp"
)

target_include_directories(ausaxs_api PUBLIC "${CMAKE_SOURCE_DIR}/include/api")
target_link_libraries(ausaxs_api ausaxs_core ausaxs_rigidbody ausaxs_em)
if (NOT EMSCRIPTEN)
	if(USE_SYSTEM_CLI11)
		find_package(CLI11 REQUIRED)
	else()
		FetchContent_Declare(
			CLI11
			GIT_REPOSITORY https://github.com/CLIUtils/CLI11
			GIT_TAG cbd58a3696887b34c70949aef21a71735a0c2ad5 # v2.7.2
			GIT_PROGRESS TRUE
		)
		FetchContent_MakeAvailable(CLI11)
	endif()
	target_sources(ausaxs_api PRIVATE
		"cli/cli_ausaxs.cpp"
		"cli/cli_debye.cpp"
		"cli/cli_em_fitter.cpp"
		"cli/cli_gpu.cpp"
		"cli/cli_rigidbody.cpp"
		"cli/cli_saxs_fitter.cpp"
		"cli/detail/CLICommon.cpp"
	)
	target_link_libraries(ausaxs_api CLI11::CLI11)
endif()
set_target_properties(ausaxs_api PROPERTIES
	POSITION_INDEPENDENT_CODE 1
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)
