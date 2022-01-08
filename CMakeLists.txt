cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Ofast -mfpmath=sse -march=native -lboost_program_options -Wall -Wextra") #-flto
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

project(SAXS)

############################################
##           Link libraries               ##
############################################
set(LIBRARIES "")
set(sources "")

# ROOT
find_package(ROOT REQUIRED)
link_directories(${ROOT_LIBRARY_DIR})

# Boost
find_package(Boost COMPONENTS program_options REQUIRED)

include_directories(${ROOT_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories("${CMAKE_SOURCE_DIR}/source")
include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories("${CMAKE_SOURCE_DIR}/test")
include_directories("${CMAKE_SOURCE_DIR}/external")

list(APPEND LIBRARIES ${ROOT_LIBRARIES})
link_libraries(${LIBRARIES})

############################################
##           Include sources              ##
############################################
include(source/io/CMakeLists.txt)
include(source/math/CMakeLists.txt)
include(source/data/CMakeLists.txt)
include(source/hydrate/CMakeLists.txt)
include(source/CMakeLists.txt)
include(source/gui/CMakeLists.txt)

############################################
##             Unit Testing               ##
############################################
message(STATUS "Compiling unit test")
file(GLOB_RECURSE TEST_SRC
  "test/*.h"
  "test/*.cpp")
add_executable(test ${TEST_SRC} ${sources})
target_link_libraries(test ${LIBRARIES})
############################################

add_subdirectory(source/scripts)
