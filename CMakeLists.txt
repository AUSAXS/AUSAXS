set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_minimum_required(VERSION 3.15)
project(AUSAXS VERSION 1.2.12)
option(GUI "Enable GUI executables" OFF)
option(DLIB "Use the dlib minimizers" ON)
option(BUILD_PLOT_EXE "Compile the plotting utility as an executable for the current platform" OFF)
option(BUILD_TESTS "Build the test suite" ON)
option(GPU "Enable WebGPU support" OFF)
option(BUILD_EXECUTABLES "Build the command-line executables" ON)
option(STATIC_CURL "Link against a static CURL library on Windows" ON)
set(ARCH "auto" CACHE STRING "Target architecture (auto, native, or a specific -march value)")
option(USE_SYSTEM_DLIB "Use system provided dlib" OFF)
option(USE_SYSTEM_GCEM "Use system provided GCEM" OFF)
option(USE_SYSTEM_BACKWARD "Use system provided backward" OFF)
option(USE_SYSTEM_CLI11 "Use system provided CLI11" OFF)
option(USE_SYSTEM_THREADPOOL "Use system provided thread-pool" OFF)
option(USE_SYSTEM_CATCH "Use system provided Catch" OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
include(cmake/setup_compile_commands.cmake)
setup_compile_commands()

############################################
##      Temporary External Library Fixes ##
############################################
# TODO: Remove once dlib fixes the missing template argument issue
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "17.0")
    add_compile_options(-Wno-error=missing-template-arg-list-after-template-kw)
endif()

############################################
##            Dependencies                ##
############################################
include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
get_property(
    compile_options
    DIRECTORY
    PROPERTY COMPILE_OPTIONS
)

set_property(
    DIRECTORY
    APPEND
    PROPERTY COMPILE_OPTIONS -w
)

if (DLIB)
	if(USE_SYSTEM_DLIB)
		find_package(dlib REQUIRED)
	else()
		FetchContent_Declare(
			dlib
			GIT_REPOSITORY https://github.com/davisking/dlib
			GIT_TAG v19.24.4
			GIT_PROGRESS TRUE
		)

		# tell dlib not to link with various unnecessary libraries
		set(DLIB_NO_GUI_SUPPORT TRUE)
		set(DLIB_JPEG_SUPPORT OFF)
		set(DLIB_LINK_WITH_SQLITE3 OFF)
		set(DLIB_USE_BLAS OFF)
		set(DLIB_USE_LAPACK OFF)
		set(DLIB_USE_CUDA OFF)
		set(DLIB_PNG_SUPPORT OFF)
		set(DLIB_GIF_SUPPORT OFF)
		set(DLIB_WEBP_SUPPORT OFF)
		set(DLIB_JXL_SUPPORT OFF)
		set(DLIB_USE_FFTW OFF)
		set(DLIB_USE_MKL_FFT OFF)
		set(DLIB_USE_FFMPEG OFF)
		set(CMAKE_CXX_STANDARD 17) # dlib must be compiled with C++17
		FetchContent_MakeAvailable(dlib)
		set(CMAKE_CXX_STANDARD 20) # continue with C++20
	endif()
	add_compile_definitions("DLIB_AVAILABLE")
endif()

if (GPU)
	FetchContent_Declare(
		webgpu
		GIT_REPOSITORY https://github.com/eliemichel/WebGPU-distribution
		GIT_TAG v0.3.0-beta
		GIT_PROGRESS TRUE
	)
	set(WEBGPU_BACKEND "DAWN")
	FetchContent_MakeAvailable(webgpu)
endif()

if(USE_SYSTEM_BACKWARD)
	find_package(Backward REQUIRED)
	find_path(BACKWARD_INCLUDE_DIR backward.hpp PATH_SUFFIXES backward)
	if (NOT BACKWARD_INCLUDE_DIR)
		message(FATAL_ERROR "backward-cpp headers not found")
	endif()
else()
	FetchContent_Declare(
		backward
		GIT_REPOSITORY https://github.com/bombela/backward-cpp
	)
	FetchContent_MakeAvailable(Backward)
	set(BACKWARD_INCLUDE_DIR ${backward_SOURCE_DIR})
endif()

if(USE_SYSTEM_THREADPOOL)
	find_path(BS_THREAD_POOL_INCLUDE_DIR BS_thread_pool.hpp)
	if (NOT BS_THREAD_POOL_INCLUDE_DIR)
		message(FATAL_ERROR "BS.thread_pool headers not found")
	endif()
else()
	FetchContent_Declare(
		thread_pool
		GIT_REPOSITORY https://github.com/bshoshany/thread-pool
	)
	FetchContent_MakeAvailable(thread_pool)
	set(BS_THREAD_POOL_INCLUDE_DIR ${thread_pool_SOURCE_DIR}/include)
endif()

if(USE_SYSTEM_GCEM)
	find_package(gcem REQUIRED)
	find_path(GCEM_INCLUDE_DIR gcem.hpp)
	if (NOT GCEM_INCLUDE_DIR)
		message(FATAL_ERROR "GCEM headers not found")
	endif()
else()
	FetchContent_Declare(
		gcem
		GIT_REPOSITORY https://github.com/klytje/gcem
	)
	FetchContent_MakeAvailable(GCEM)
	set(GCEM_INCLUDE_DIR ${gcem_SOURCE_DIR}/include)
endif()

set_property(
    DIRECTORY
    PROPERTY COMPILE_OPTIONS ${compile_options}
)
unset(compile_options)

include_directories(${BS_THREAD_POOL_INCLUDE_DIR} ${GCEM_INCLUDE_DIR} ${BACKWARD_INCLUDE_DIR})

############################################
##           Find and link CURL           ##
############################################
if (WIN32)
	if (STATIC_CURL)
		add_compile_definitions("CURL_STATICLIB")
	endif()
	find_package(CURL REQUIRED)
	set(LIBS -static CURL::libcurl)
	include_directories(${CURL_INCLUDE_DIRS}) # find_package apparently does not always include the headers
elseif (APPLE)
	find_package(CURL REQUIRED)
	set(LIBS CURL::libcurl)
	link_libraries("$<$<CONFIG:Debug>:-ldwarf>")
elseif(UNIX)
	find_package(CURL REQUIRED)
	set(LIBS CURL::libcurl -static-libgcc -static-libstdc++)
	link_libraries("$<$<CONFIG:Debug>:-ldwarf>")
endif()
############################################
##                Doxygen                 ##
############################################
find_package(Doxygen)
if (DOXYGEN_FOUND)
	set(sim3a_Doxygen "${CMAKE_BINARY_DIR}/saxs.dox")
	configure_file(${CMAKE_SOURCE_DIR}/scripts/ausaxs.dox.in ${sim3a_Doxygen} @ONLY)
	add_custom_target(
		doc
		${DOXYGEN_EXECUTABLE} ${sim3a_Doxygen}
    		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)
endif()


############################################
##           Build library                ##
############################################
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
	include(CheckIPOSupported)
	check_ipo_supported(RESULT LTO_SUPPORTED)
else()
	set(LTO_SUPPORTED OFF)
endif()

add_subdirectory(source)
if (BUILD_TESTS)
	add_subdirectory(tests)
endif()
if (BUILD_EXECUTABLES)
	add_subdirectory(executable)
endif()
