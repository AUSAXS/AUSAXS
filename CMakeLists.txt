cmake_minimum_required(VERSION 3.14)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE DEBUG)

project(SAXS VERSION 0.0)
add_compile_definitions("$<$<CONFIG:DEBUG>:DEBUG=1;SAFE_MATH=1>")#;_GLIBCXX_DEBUG>")
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
FetchContent_MakeAvailable(CLI11 backward BSThreadPool)

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
file(GLOB_RECURSE HEADER_FILES "${CMAKE_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE SOURCE_FILES "${CMAKE_SOURCE_DIR}/source/*.cpp")

add_library(objlib OBJECT ${SOURCE_FILES})
set_target_properties(objlib PROPERTIES POSITION_INDEPENDENT_CODE 1)

add_library(saxs_static STATIC $<TARGET_OBJECTS:objlib>)
add_library(saxs SHARED $<TARGET_OBJECTS:objlib>)
set_target_properties(saxs_static PROPERTIES OUTPUT_NAME saxs_static)
set_target_properties(saxs PROPERTIES OUTPUT_NAME saxs)

target_include_directories(objlib PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
target_include_directories(saxs PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
target_include_directories(saxs_static PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
target_link_libraries(objlib CLI11::CLI11 dlib::dlib ${CURL})
target_link_libraries(saxs CLI11::CLI11 dlib::dlib ${CURL})
target_link_libraries(saxs_static CLI11::CLI11 dlib::dlib ${CURL})
link_libraries(saxs_static)

############################################
##             Unit Testing               ##
############################################
add_subdirectory(test)

############################################
##          Other executables             ##
############################################
add_subdirectory(executable)
