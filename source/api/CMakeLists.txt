if(USE_SYSTEM_CLI11)
	find_package(CLI11 REQUIRED)
else()
	FetchContent_Declare(
		CLI11
		GIT_REPOSITORY https://github.com/CLIUtils/CLI11
		GIT_TAG 792d89286788acac125e0487f8dbde88035f7422
		GIT_PROGRESS TRUE
	)
	FetchContent_MakeAvailable(CLI11)
endif()

add_library(ausaxs_api OBJECT 
    "ObjectStorage.cpp"

    "api_helper.cpp"
    "api_sasview.cpp"

    "cli/cli_em_fitter.cpp"
    "cli/cli_rigidbody.cpp"
    "cli/cli_saxs_fitter.cpp"

    "pyausaxs/api_data.cpp"
    "pyausaxs/api_fit.cpp"
    "pyausaxs/api_io.cpp"
    "pyausaxs/api_iterative_fit.cpp"
    "pyausaxs/api_molecule.cpp"
    "pyausaxs/api_pdb.cpp"
    "pyausaxs/api_settings.cpp"
)

target_include_directories(ausaxs_api PUBLIC "${CMAKE_SOURCE_DIR}/include/api")
target_link_libraries(ausaxs_api ausaxs_core ausaxs_rigidbody ausaxs_em CLI11::CLI11)
set_target_properties(ausaxs_api PROPERTIES POSITION_INDEPENDENT_CODE 1)
