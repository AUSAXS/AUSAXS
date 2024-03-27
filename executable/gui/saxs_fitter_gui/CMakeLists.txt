project(SAXSFitterGUI)

set(ELEMENTS_APP_PROJECT "saxs_fitter_gui")
set(ELEMENTS_APP_TITLE "SAXS fitter GUI")
set(ELEMENTS_APP_COPYRIGHT "Copyright (c) 2024 Kristian Lytje")
set(ELEMENTS_APP_ID ".com.ausaxs.saxs-fitter")
set(ELEMENTS_APP_VERSION "1.0")
set(ELEMENTS_APP_SOURCES ${sources} ${CMAKE_CURRENT_SOURCE_DIR}/saxs_fitter_gui.cpp)
include(ElementsConfigApp)

target_link_libraries("${ELEMENTS_APP_PROJECT}" PRIVATE ausaxs_static)
set_target_properties("${ELEMENTS_APP_PROJECT}" PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
if (BUILD_PLOT_EXE)
	add_dependencies("${ELEMENTS_APP_PROJECT}" plotter)
endif()

add_custom_command(TARGET saxs_fitter_gui POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/scripts" "${CMAKE_BINARY_DIR}/executable/gui/em_fitter_gui/scripts"
)