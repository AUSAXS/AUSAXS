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
target_link_libraries("${ELEMENTS_APP_PROJECT}" PRIVATE ausaxs_core ausaxs_math ausaxs_em ausaxs_gui nfd)
target_include_directories("${ELEMENTS_APP_PROJECT}" PRIVATE "${CMAKE_SOURCE_DIR}/include/gui")
add_dependencies("${ELEMENTS_APP_PROJECT}" plotter)
