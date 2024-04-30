FetchContent_Declare(
	CLI11
	GIT_REPOSITORY https://github.com/CLIUtils/CLI11
	GIT_TAG 792d89286788acac125e0487f8dbde88035f7422
	GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(CLI11)

file(GLOB EXEC_SOURCE "*.cpp")
foreach(exec ${EXEC_SOURCE})
	cmake_path(GET exec STEM name)
	add_executable("${name}" ${exec})
	set_target_properties("${name}" PROPERTIES 
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
	)
	target_link_libraries("${name}" CLI11::CLI11 ausaxs_static "$<$<CONFIG:Debug>:Backward::Interface>")
endforeach()

if (GUI)
	add_subdirectory(gui)
endif()
