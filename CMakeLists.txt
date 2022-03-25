cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Ofast -mfpmath=sse -march=native -lboost_program_options -Wall -Wextra -g") #-flto
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

project(SAXS)


############################################
##           Link libraries               ##
############################################
set(LIBRARIES)

# ROOT
find_package(ROOT REQUIRED)
link_directories(${ROOT_LIBRARY_DIR})

# Boost
find_package(Boost COMPONENTS program_options REQUIRED)

include_directories(${ROOT_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

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
file(GLOB_RECURSE HEADERS "include/*.h")
file(GLOB_RECURSE SOURCES "source/*.cpp")
add_library(saxs ${HEADERS} ${SOURCES})
target_include_directories (saxs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/source")


############################################
##             Unit Testing               ##
############################################
message(STATUS "Compiling unit test")
file(GLOB_RECURSE TEST_SRC
  "test/*.h"
  "test/*.cpp")
add_executable(test ${TEST_SRC})
target_link_libraries(test saxs)
target_include_directories(test PRIVATE "${CMAKE_SOURCE_DIR}/test" "${CMAKE_CURRENT_SOURCE_DIR}/external" "${CMAKE_CURRENT_SOURCE_DIR}external/CLI11.hpp")


############################################
##          Other executables             ##
############################################
add_subdirectory(executable)
