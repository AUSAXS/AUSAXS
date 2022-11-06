cmake_minimum_required(VERSION 3.14)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Ofast -mfpmath=sse -march=native -Wall -Wpedantic -Wextra -g -fno-finite-math-only")# --coverage -g")# -fsanitize=address") #-flto
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CXX_STANDARD C++17)

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

FetchContent_Declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11
  GIT_TAG b9be5b9444772324459989177108a6a65b8b2769
  GIT_PROGRESS TRUE
)

FetchContent_Declare(
  dlib
  GIT_REPOSITORY https://github.com/davisking/dlib
  GIT_TAG v19.18
  GIT_PROGRESS TRUE
)
set(DLIB_ENABLE_ASSERTS false)

FetchContent_MakeAvailable(Catch2 CLI11 dlib)


############################################
##           Find and link CURL           ##
############################################
# when building for windows we want a single executable, so we have to include a bunch of static libraries to make curl work
if (WIN32)
	add_definitions(-DCURL_STATICLIB)
	set(OPENSSL_USE_STATIC_LIBS TRUE)
	set(CURL_LIBRARY "/usr/x86_64-w64-mingw32/usr/local/lib/libcurl.a")
	find_package(CURL REQUIRED)
	find_package(OpenSSL COMPONENTS Crypto SSL REQUIRED)
	set(CURL -static gcc stdc++ winpthread CURL::libcurl OpenSSL::SSL OpenSSL::Crypto ws2_32 libwldap32.a crypt32)

# for linux we assume the users can figure out how to install missing packages by themselves
elseif (UNIX)
	find_package(CURL REQUIRED)
	find_package(OpenSSL COMPONENTS Crypto SSL REQUIRED)
	set(CURL CURL::libcurl OpenSSL::SSL OpenSSL::Crypto)
endif()


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
##           Build library                ##
############################################
file(GLOB_RECURSE HEADERS "include/*.h")
file(GLOB_RECURSE SOURCES "source/*.cpp")
add_library(saxs ${HEADERS} ${SOURCES})
target_include_directories (saxs PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/source")
target_link_libraries(saxs CLI11::CLI11 dlib::dlib ${CURL})


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
