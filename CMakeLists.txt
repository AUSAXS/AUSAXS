cmake_minimum_required(VERSION 3.14)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Ofast -mfpmath=sse -march=native -lboost_program_options -Wall -Wextra -g -fno-finite-math-only")# --coverage -g")# -fsanitize=address") #-flto
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

project(SAXS VERSION 0.0)


############################################
##            Dependencies                ##
############################################
include(FetchContent)
#set(FETCHCONTENT_UPDATES_DISCONNECT ON)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG f4af9f69265d009a457aa99d1075cfba78652a66
  GIT_PROGRESS TRUE
)

include(FetchContent)
FetchContent_Declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11
  GIT_TAG b9be5b9444772324459989177108a6a65b8b2769
  GIT_PROGRESS TRUE
)

FetchContent_MakeAvailable(Catch2 CLI11)


############################################
##           Link libraries               ##
############################################
set(LIBRARIES)

# ROOT
find_package(ROOT REQUIRED)
link_directories(${ROOT_LIBRARY_DIR})
include_directories(${ROOT_INCLUDE_DIR})

list(APPEND LIBRARIES ${ROOT_LIBRARIES})
link_libraries(${LIBRARIES})
link_libraries(CLI11::CLI11)

find_package(CURL REQUIRED)
link_libraries(CURL::libcurl)


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
file(GLOB_RECURSE TEST_SRC
  "test/*.h"
  "test/*.cpp")
add_executable(test ${TEST_SRC})
target_link_libraries(test PRIVATE saxs Catch2::Catch2WithMain)


############################################
##          Other executables             ##
############################################
add_subdirectory(executable)
