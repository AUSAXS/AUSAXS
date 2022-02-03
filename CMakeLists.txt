cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Ofast -mfpmath=sse -march=native -lboost_program_options -Wall -Wextra -g") #-flto
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra -g") #-flto
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

project(SAXS)

############################################
##           Link libraries               ##
############################################
set(LIBRARIES)
set(sources "${CMAKE_SOURCE_DIR}/external/CLI11.hpp")

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
##                Doxygen                 ##
############################################
find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(sim3a_Doxygen "${CMAKE_CURRENT_BINARY_DIR}/saxs.dox")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/saxs.dox.in ${sim3a_Doxygen} @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${sim3a_Doxygen}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

############################################
##           Include sources              ##
############################################
include(source/CMakeLists.txt)
include(source/io/CMakeLists.txt)
include(source/math/CMakeLists.txt)
include(source/data/CMakeLists.txt)
include(source/rigidbody/CMakeLists.txt)
include(source/hydrate/CMakeLists.txt)
include(source/fitter/CMakeLists.txt)
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
