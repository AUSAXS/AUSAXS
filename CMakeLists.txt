set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_minimum_required(VERSION 3.15)
project(AUSAXS VERSION 1.2.0)
option(GUI "Enable GUI executables" OFF)
option(DLIB "Download and use the dlib minimizers" ON)
option(BUILD_PLOT_EXE "Compile the plotting utility as an executable for the current platform" OFF)
option(CONSTEXPR_TABLES "Generate lookup tables at compile-time" OFF)
set(ARCH "native" CACHE STRING "Target architecture. Default: native")
option(USE_SYSTEM_GCEM "Use system provided GCEM" OFF)
option(USE_SYSTEM_BACKWARD "Use system provided backward" OFF)
option(USE_SYSTEM_CLI11 "Use system provided CLI11" OFF)
option(USE_SYSTEM_THREADPOOL "Use system provided thread-pool" OFF)
option(USE_SYSTEM_CATCH "Use system provided Catch" OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
add_compile_definitions("CONSTEXPR_TABLES=${CONSTEXPR_TABLES};$<$<CONFIG:DEBUG>:DEBUG=1;SAFE_MATH=1>")
if (WIN32)
	if (MSVC)
		add_compile_definitions("NOMINMAX;BUILD_EXPORT_DLL")
		add_compile_options(
			/fp:fast /constexpr:steps10000000000 /Zm500 
			/wd4267 # disable size_t --> int, unsigned int conversions
			/wd4244 # disable double --> float,int conversions
			"$<$<STREQUAL:${ARCH},x86-64>:/arch:AVX>"
			"$<$<STREQUAL:${ARCH},arm64>:/arch:armv8.0>"
		)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
	elseif (MINGW)
		add_compile_options(
			-O3 -ffast-math -fno-finite-math-only -pipe -mavx
			"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},Clang>:-fconstexpr-steps=1000000000>"
			"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:-fconstexpr-ops-limit=10000000000>"
			"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra;-march=native>"
			"$<$<AND:$<CONFIG:RELEASE>,$<STREQUAL:${ARCH},x86-64>>:-march=x86-64-v3>"
			"$<$<AND:$<CONFIG:RELEASE>,$<STREQUAL:${ARCH},native>>:-march=native>"
		)
	endif()
elseif (APPLE)
	add_compile_options(
		-O3 -ffast-math -fno-finite-math-only -pipe 
		"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},Clang>:-fconstexpr-steps=1000000000>"
		"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:-fconstexpr-ops-limit=10000000000>"
		"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra>"
		"$<$<AND:$<CONFIG:RELEASE>,$<STREQUAL:${ARCH},x86-64>>:-march=x86-64>"
		"$<$<AND:$<CONFIG:RELEASE>,$<STREQUAL:${ARCH},arm64>>:-march=armv8-a>"
	)
elseif (UNIX)
#	add_compile_definitions("$<$<CONFIG:DEBUG>:_GLIBCXX_DEBUG>")
	add_compile_options(
		-Ofast -fno-finite-math-only -pipe -mavx
		"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},Clang>:-fconstexpr-steps=1000000000>"
		"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:-fconstexpr-ops-limit=10000000000>"
		"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra;-march=native>"
		"$<$<AND:$<CONFIG:RELEASE>,$<STREQUAL:${ARCH},x86-64>>:-march=x86-64-v3>"
		"$<$<AND:$<CONFIG:RELEASE>,$<STREQUAL:${ARCH},native>>:-march=native>"
	)
endif()

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
	add_compile_definitions("DLIB_AVAILABLE")
endif()

if(USE_SYSTEM_BACKWARD)
	find_package(Backward REQUIRED)
else()
	FetchContent_Declare(
		backward
		GIT_REPOSITORY https://github.com/bombela/backward-cpp
	)
	FetchContent_MakeAvailable(Backward)
endif()

if(USE_SYSTEM_THREADPOOL)
	include(CheckIncludeFileCXX)

	check_include_file_cxx(
		BS_thread_pool.hpp
		BS_THREAD_POOL_FOUND
	)

	if (NOT BS_THREAD_POOL_FOUND)
		message(FATAL_ERROR "BS.thread_pool headers not found")
	endif()
else()
	FetchContent_Declare(
		thread_pool
		GIT_REPOSITORY https://github.com/bshoshany/thread-pool
	)
	FetchContent_MakeAvailable(thread_pool)
endif()

if(USE_SYSTEM_GCEM)
	find_package(gcem REQUIRED)
else()
	FetchContent_Declare(
		gcem
		GIT_REPOSITORY https://github.com/klytje/gcem
	)
	FetchContent_MakeAvailable(GCEM)
endif()

set_property(
    DIRECTORY
    PROPERTY COMPILE_OPTIONS ${compile_options}
)
unset(compile_options)

include_directories(${thread_pool_SOURCE_DIR}/include ${gcem_SOURCE_DIR}/include ${backward_SOURCE_DIR})

############################################
##           Find and link CURL           ##
############################################
if (WIN32)
	add_compile_definitions("CURL_STATICLIB")
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
add_subdirectory(tests)
add_subdirectory(executable)
