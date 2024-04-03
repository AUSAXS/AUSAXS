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

# plotter exe
set(BUILD_PLOT_EXE FALSE)
if (WIN32)
	find_package(Python3 COMPONENTS Interpreter)
	if (Python3_FOUND)
		execute_process(
			COMMAND pip show pyinstaller matplotlib numpy
			RESULT_VARIABLE EXIT_CODE
			OUTPUT_QUIET
		)
		if (${EXIT_CODE} EQUAL 0)
			add_custom_target(plotter)
			add_custom_command(TARGET plotter POST_BUILD
				COMMAND ${Python3_EXECUTABLE} -c "import PyInstaller.__main__; PyInstaller.__main__.run(['${CMAKE_SOURCE_DIR}/scripts/plot.py', '--onefile', '--distpath', '${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bin'])"
				VERBATIM
			)
			set(BUILD_PLOT_EXE TRUE)
		endif()
	endif()
endif()

# gui executables
add_subdirectory(saxs_fitter_gui)
add_subdirectory(em_fitter_gui)
add_subdirectory(saxs_fitter_gui)
