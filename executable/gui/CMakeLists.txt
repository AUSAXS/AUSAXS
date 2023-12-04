###############################################################################
#  Copyright (c) 2016-2020 Joel de Guzman
#
#  Distributed under the MIT License (https://opensource.org/licenses/MIT)
###############################################################################
cmake_minimum_required(VERSION 3.9.6...3.15.0)
project(Buttons LANGUAGES C CXX)

FetchContent_Declare(
	elements
	GIT_REPOSITORY https://github.com/cycfi/elements.git
)

set(CMAKE_CXX_STANDARD 20) # continue with C++20
FetchContent_MakeAvailable(CLI11 BSThreadPool elements)

set(ELEMENTS_ROOT "${elements_SOURCE_DIR}")
get_filename_component(ELEMENTS_ROOT "${ELEMENTS_ROOT}" ABSOLUTE)
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${ELEMENTS_ROOT}/cmake")

# If we are building outside the project, you need to set ELEMENTS_ROOT:
if (NOT ELEMENTS_BUILD_EXAMPLES)
#   include(ElementsConfigCommon)
   set(ELEMENTS_BUILD_EXAMPLES OFF)
   add_subdirectory(${ELEMENTS_ROOT} elements)
endif()

set(ELEMENTS_APP_PROJECT "gui")
set(ELEMENTS_APP_TITLE "gui")
set(ELEMENTS_APP_COPYRIGHT "Copyright (c) 2016-2020 Joel de Guzman")
set(ELEMENTS_APP_ID "gui")
set(ELEMENTS_APP_VERSION "1.0")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(ELEMENTS_APP_SOURCES ${sources} ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/range_slider.cpp)
#set(ELEMENTS_APP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/gui/main.cpp)

include(ElementsConfigApp)
