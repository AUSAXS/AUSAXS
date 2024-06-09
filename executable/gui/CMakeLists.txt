FetchContent_Declare(
	elements
	GIT_REPOSITORY https://github.com/cycfi/elements.git
)
FetchContent_Declare(
	nfd
	GIT_REPOSITORY https://github.com/btzy/nativefiledialog-extended.git
)
FetchContent_MakeAvailable(elements CLI11 nfd)

link_libraries(nfd)
include_directories(include)

# set the path to the elements root & module path
set(ELEMENTS_ROOT "${elements_SOURCE_DIR}")
set(ELEMENTS_BUILD_EXAMPLES OFF)
get_filename_component(ELEMENTS_ROOT "${ELEMENTS_ROOT}" ABSOLUTE)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${ELEMENTS_ROOT}/cmake")

# target for ensuring plotting is possible for the gui executables
add_custom_target(plotter)
if (${BUILD_PLOT_EXE})
	find_package(Python3 COMPONENTS Interpreter)
	if (Python3_FOUND)
		execute_process(
			COMMAND pip show pyinstaller matplotlib numpy
			RESULT_VARIABLE EXIT_CODE
			OUTPUT_QUIET
		)
		if (${EXIT_CODE} EQUAL 0)
		  add_custom_command(TARGET plotter POST_BUILD
				COMMAND ${Python3_EXECUTABLE} -c "import PyInstaller.__main__; PyInstaller.__main__.run(['${CMAKE_SOURCE_DIR}/scripts/plot.py', '--onefile', '--distpath', '${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bin'])"
				VERBATIM
			)
		elseif(${GUI})
			message(FATAL_ERROR "GUI executables requested, but Python is not installed. The plotting script is unusable. Compile using GUI=FALSE to avoid this error.")
		else()
			message(WARNING "The plotting utility cannot be built. The GUI executables will only be able to use the raw Python scripts.")
		endif()
	endif()
elseif(${GUI})
	find_package(Python3 COMPONENTS Interpreter)
	add_custom_command(TARGET plotter POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/scripts" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/scripts"
	)
endif()

# gui executables
add_subdirectory(saxs_fitter_gui)
add_subdirectory(em_fitter_gui)
