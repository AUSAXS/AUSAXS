cmake_minimum_required(VERSION 3.14)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE RELEASE)

project(SAXS VERSION 0.0)
add_compile_definitions("$<$<CONFIG:DEBUG>:DEBUG=1;SAFE_MATH=1>")
add_compile_options(
	-Ofast -fno-finite-math-only# --coverage
	"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra>"
	"$<$<CONFIG:RELEASE>:-flto;-march=native>"
)


############################################
##            Dependencies                ##
############################################
include(FetchContent)
#set(FETCHCONTENT_UPDATES_DISCONNECT ON)
FetchContent_Declare(
	dlib
	GIT_REPOSITORY https://github.com/davisking/dlib
	GIT_TAG v19.18
	GIT_PROGRESS TRUE
)
set(DLIB_NO_GUI_SUPPORT TRUE)
FetchContent_MakeAvailable(dlib)

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
	backward
	GIT_REPOSITORY https://github.com/bombela/backward-cpp
	GIT_TAG v1.6
)

FetchContent_Declare(
	BSThreadPool
	GIT_REPOSITORY https://github.com/darren-harton/thread-pool
)
set(CMAKE_CXX_STANDARD 20)
FetchContent_MakeAvailable(Catch2 CLI11 backward BSThreadPool)

include_directories(${CMAKE_BINARY_DIR}/_deps/bsthreadpool-src)

# we cannot debug without debug symbols
set(ENABLE_DEBUG $<$<CONFIG:DEBUG>: ${BACKWARD_ENABLE}>)


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
#	if (NOT BUILD_SHARED_LIBS)
#		set(OPENSSL_USE_STATIC_LIBS TRUE)
#		set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#	endif()
	find_package(CURL REQUIRED)
	find_package(OpenSSL COMPONENTS Crypto SSL REQUIRED)
#	link_libraries(-static -static-libgcc -static-libstdc++)
#	set(CURL CURL::libcurl OpenSSL::SSL OpenSSL::Crypto)
	set(CURL -static-libgcc -static-libstdc++ CURL::libcurl OpenSSL::SSL OpenSSL::Crypto)
	if (CMAKE_BUILD_TYPE MATCHES DEBUG)
		link_libraries(-lbfd -ldl -lgcov)
	endif()
endif()


############################################
##                Doxygen                 ##
############################################
find_package(Doxygen)
if(DOXYGEN_FOUND)
	set(sim3a_Doxygen "${CMAKE_BINARY_DIR}/saxs.dox")
	configure_file(${CMAKE_SOURCE_DIR}/saxs.dox.in ${sim3a_Doxygen} @ONLY)
	add_custom_target(
		doc
		${DOXYGEN_EXECUTABLE} ${sim3a_Doxygen}
    		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)
endif(DOXYGEN_FOUND)


############################################
##           Build library                ##
############################################
file(GLOB_RECURSE HEADERS "${CMAKE_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/source/*.cpp")
add_library(saxs ${HEADERS} ${SOURCES})
target_include_directories(saxs PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
target_link_libraries(saxs CLI11::CLI11 dlib::dlib ${CURL})
link_libraries(saxs)

############################################
##             Unit Testing               ##
############################################
file(GLOB_RECURSE 
	TEST_SRC
	"test/*.h"
	"test/*.cpp"
)

# make the tests available through cmake
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)

add_custom_target(tests)
foreach(test ${TEST_SRC})
	cmake_path(GET test STEM testname)
	set(name "test_${testname}")
	add_executable("${name}" ${test} ${ENABLE_DEBUG})
	target_link_libraries("${name}" Catch2::Catch2WithMain)
	set_target_properties("${name}" PROPERTIES 
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
	)
	if (CMAKE_BUILD_TYPE MATCHES DEBUG)
		add_backward("${name}")
	endif()	
	add_dependencies(tests "${name}")
	catch_discover_tests("${name}" TEST_SPEC "~[slow] ~[broken]" WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/..")
endforeach()

############################################
##          Other executables             ##
############################################
add_subdirectory(executable)
