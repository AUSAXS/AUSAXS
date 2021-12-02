cmake_minimum_required(VERSION 2.8)
project(Calibrate)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -flto -lboost_program_options")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

#----------------------------------------------------------------------------
set(LIBRARIES "")
set(sources "")

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
#----------------------------------------------------------------------------
include(source/math/CMakeLists.txt)
include(source/data/CMakeLists.txt)
include(source/hydrate/CMakeLists.txt)
include(source/CMakeLists.txt)

add_subdirectory(source/tests)
add_subdirectory(source/scripts)
