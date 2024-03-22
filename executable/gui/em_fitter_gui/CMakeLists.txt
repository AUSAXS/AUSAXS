project(EMGUI)

set(ELEMENTS_APP_PROJECT "em_fitter_gui")
set(ELEMENTS_APP_TITLE "EM fitter GUI")
set(ELEMENTS_APP_COPYRIGHT "Copyright (c) 2024 Kristian Lytje")
set(ELEMENTS_APP_ID ".com.ausaxs.em-fitter")
set(ELEMENTS_APP_VERSION "1.0")
set(ELEMENTS_APP_SOURCES ${sources} ${CMAKE_CURRENT_SOURCE_DIR}/em_fitter_gui.cpp)
include(ElementsConfigApp)

add_custom_command(TARGET em_fitter_gui POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/scripts" "${CMAKE_BINARY_DIR}/executable/gui/em_fitter_gui/scripts"
)
