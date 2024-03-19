FetchContent_Declare(
	elements
	GIT_REPOSITORY https://github.com/cycfi/elements.git
	GIT_TAG develop
)

FetchContent_Declare(
	nfd
	GIT_REPOSITORY https://github.com/btzy/nativefiledialog-extended.git
)

FetchContent_MakeAvailable(CLI11 nfd)
link_libraries(nfd)
include_directories(include)

project(GUI)
set(ELEMENTS_ROOT "${elements_SOURCE_DIR}")
set(ELEMENTS_BUILD_EXAMPLES OFF)

get_filename_component(ELEMENTS_ROOT "${ELEMENTS_ROOT}" ABSOLUTE)
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${ELEMENTS_ROOT}/cmake")

add_subdirectory(${ELEMENTS_ROOT} elements)
add_subdirectory(saxs_fitter_gui)
add_subdirectory(em_fitter_gui)
