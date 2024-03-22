project(EMGUI)

include(ElementsConfigCommon)
set(ELEMENTS_APP_PROJECT "em_fitter_gui")
set(ELEMENTS_APP_TITLE "EM fitter GUI")
set(ELEMENTS_APP_COPYRIGHT "Copyright (c) 2024 Kristian Lytje")
set(ELEMENTS_APP_ID ".com.ausaxs.em-fitter")
set(ELEMENTS_APP_VERSION "1.0")
set(ELEMENTS_APP_SOURCES ${sources} ${CMAKE_CURRENT_SOURCE_DIR}/em_fitter_gui.cpp)
include(ElementsConfigApp)

#add_library(ausaxs_imp SHARED IMPORTED)
#set_target_properties(ausaxs_imp PROPERTIES
#	IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/Release/ausaxs.dll"
#	IMPORTED_IMPLIB "${CMAKE_BINARY_DIR}/Release/ausaxs.lib"
#	INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include"
#)
target_link_libraries("${ELEMENTS_APP_PROJECT}" PRIVATE ausaxs_static)
target_include_directories("${ELEMENTS_APP_PROJECT}" PUBLIC "${CMAKE_SOURCE_DIR}/include")
set_target_properties("${ELEMENTS_APP_PROJECT}" PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
if (BUILD_PLOT_EXE)
	add_dependencies("${ELEMENTS_APP_PROJECT}" plotter)
endif()

add_custom_command(TARGET em_fitter_gui POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/scripts" "${CMAKE_BINARY_DIR}/executable/gui/em_fitter_gui/scripts"
<<<<<<< HEAD
)
=======
)
>>>>>>> 30c39a5d (cmake fixes)
