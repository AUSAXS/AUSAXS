cmake_minimum_required(VERSION 2.8)
project(Calibrate)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -lboost_program_options")
SET (CMAKE_CXX_COMPILER "/usr/bin/g++")

#----------------------------------------------------------------------------
# Dependencies
#

# Variable for holding dependant libraries
set(LIBRARIES "")

# Include ROOT
find_package(ROOT REQUIRED)
link_directories(${ROOT_LIBRARY_DIR})

# Include Boost
find_package(Boost COMPONENTS program_options REQUIRED)

include_directories(${ROOT_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories("${CMAKE_SOURCE_DIR}/source")

list(APPEND LIBRARIES ${ROOT_LIBRARIES})

link_libraries(${LIBRARIES})

add_executable(main main.cpp)

add_subdirectory(source/tests)
add_subdirectory(source/scripts)
