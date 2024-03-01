###############################################################################
#  Copyright (c) 2016-2020 Joel de Guzman
#
#  Distributed under the MIT License (https://opensource.org/licenses/MIT)
###############################################################################
cmake_minimum_required(VERSION 3.9.6...3.15.0)
project(EMGUI LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
FetchContent_Declare(
	elements
	GIT_REPOSITORY https://github.com/cycfi/elements.git
	GIT_TAG develop
)

FetchContent_Declare(
	nfd
	GIT_REPOSITORY https://github.com/btzy/nativefiledialog-extended.git
)
FetchContent_MakeAvailable(elements CLI11 nfd)
link_libraries(nfd)
include_directories(resources)

set(ELEMENTS_ROOT "${elements_SOURCE_DIR}")
get_filename_component(ELEMENTS_ROOT "${ELEMENTS_ROOT}" ABSOLUTE)
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${ELEMENTS_ROOT}/cmake")

if (NOT ELEMENTS_BUILD_EXAMPLES)
   set(ELEMENTS_BUILD_EXAMPLES OFF)
   add_subdirectory(${ELEMENTS_ROOT} elements)
endif()

set(ELEMENTS_APP_PROJECT "gui")
set(ELEMENTS_APP_TITLE "gui")
set(ELEMENTS_APP_COPYRIGHT "Copyright (c) 2023 Kristian Lytje")
set(ELEMENTS_APP_ID "gui")
set(ELEMENTS_APP_VERSION "1.0")
set(ELEMENTS_APP_SOURCES ${sources} ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

include(ElementsConfigApp)
