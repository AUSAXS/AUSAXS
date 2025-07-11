project(SAXSFitterGUI)

set(ELEMENTS_APP_PROJECT "saxs_fitter_gui")
set(ELEMENTS_APP_TITLE "SAXS fitter GUI")
set(ELEMENTS_APP_COPYRIGHT "Copyright (c) 2024 Kristian Lytje")
set(ELEMENTS_APP_ID "com.ausaxs.saxs_fitter_gui")
set(ELEMENTS_APP_VERSION "1.0")
set(ELEMENTS_APP_SOURCES ${sources} ${CMAKE_CURRENT_SOURCE_DIR}/saxs_fitter_gui.cpp)
include(ElementsConfigApp)

target_link_libraries("${ELEMENTS_APP_PROJECT}" PRIVATE ausaxs_core ausaxs_math ausaxs_gui nfd)
target_include_directories("${ELEMENTS_APP_PROJECT}" PRIVATE "${CMAKE_SOURCE_DIR}/include/gui")
add_dependencies("${ELEMENTS_APP_PROJECT}" plotter)
