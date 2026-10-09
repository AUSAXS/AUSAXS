set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_minimum_required(VERSION 3.15)
project(AUSAXS VERSION 1.3.0)
option(POCKETFFT "Use pocketfft for the lattice transform of the grid excluded volume" ON)
option(BUILD_TESTS "Build the test suite" ON)
option(BUILD_EXECUTABLES "Build the command-line executables" ON)
option(STATIC_CURL "Link against a static CURL library on Windows" ON)
set(ARCH "auto" CACHE STRING "Target architecture (auto, native, or a specific -march value)")
option(USE_SYSTEM_GCEM "Use system provided GCEM" OFF)
option(USE_SYSTEM_CLI11 "Use system provided CLI11" OFF)
option(USE_SYSTEM_THREADPOOL "Use system provided thread-pool" OFF)
option(USE_SYSTEM_POCKETFFT "Use system provided pocketfft" OFF)
option(USE_SYSTEM_CATCH "Use system provided Catch" OFF)
option(GPU_DOWNLOADER "Build the downloader for prebuilt GPU backends (ausaxs gpu install). Needs miniz." ON)
option(USE_SYSTEM_MINIZ "Use system provided miniz" OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# A GPU backend is opened by name at run time rather than linked, so it is not on any target's link
# line and nothing would otherwise teach the binaries in the build tree where to find it.
list(APPEND CMAKE_BUILD_RPATH "${CMAKE_BINARY_DIR}/lib")
include(cmake/setup_compile_commands.cmake)
setup_compile_commands()

############################################
##            Dependencies                ##
############################################
include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
if (POLICY CMP0135)
	cmake_policy(SET CMP0135 NEW) # URL downloads get the time of extraction as their timestamp
endif()
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

if (POCKETFFT)
	if(USE_SYSTEM_POCKETFFT)
		find_path(POCKETFFT_INCLUDE_DIR pocketfft_hdronly.h)
		if (NOT POCKETFFT_INCLUDE_DIR)
			message(FATAL_ERROR "pocketfft headers not found")
		endif()
	else()
		FetchContent_Declare(
			pocketfft
			GIT_REPOSITORY https://github.com/mreineck/pocketfft
			GIT_TAG c90e55b3d529f8efa40ed01a20de22405f45fc65 # the header-only C++ version lives on the "cpp" branch, which carries no tags
		)
		FetchContent_MakeAvailable(pocketfft)
		set(POCKETFFT_INCLUDE_DIR ${pocketfft_SOURCE_DIR})
	endif()
	add_compile_definitions("POCKETFFT_AVAILABLE")
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

if (GPU_DOWNLOADER)
	if (USE_SYSTEM_MINIZ)
		find_package(miniz CONFIG REQUIRED)
	else()
		# the release archive holds the amalgamated miniz.c/miniz.h, which build without miniz's own CMake configuration
		FetchContent_Declare(
			miniz
			URL https://github.com/richgel999/miniz/releases/download/3.1.2/miniz-3.1.2.zip
			URL_HASH SHA256=f0446d863f9c19926ad9483c523fdc42e42b8d4a6a431d27e09d49c79a140d9a
		)
		FetchContent_MakeAvailable(miniz)
		add_library(miniz STATIC "${miniz_SOURCE_DIR}/miniz.c")
		target_include_directories(miniz SYSTEM PUBLIC "${miniz_SOURCE_DIR}")
		set_target_properties(miniz PROPERTIES POSITION_INDEPENDENT_CODE ON)
		add_library(miniz::miniz ALIAS miniz)
	endif()
endif()

set_property(
    DIRECTORY
    PROPERTY COMPILE_OPTIONS ${compile_options}
)
unset(compile_options)
include_directories(${BS_THREAD_POOL_INCLUDE_DIR} ${GCEM_INCLUDE_DIR})
if (POCKETFFT)
	include_directories(SYSTEM ${POCKETFFT_INCLUDE_DIR})
endif()

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
elseif(UNIX)
	find_package(CURL REQUIRED)
	set(LIBS CURL::libcurl -static-libgcc -static-libstdc++)
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
string(TOUPPER "${CMAKE_BUILD_TYPE}" AUSAXS_BUILD_TYPE_UPPER)
if (NOT AUSAXS_BUILD_TYPE_UPPER STREQUAL "RELEASE")
	set(LTO_SUPPORTED OFF)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" AND BUILD_TESTS)
	set(LTO_SUPPORTED OFF)
	message(STATUS "LTO disabled: MSVC cannot combine it with the symbol export BUILD_TESTS requires")
else()
	include(CheckIPOSupported)
	check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_MESSAGE)
	if (NOT LTO_SUPPORTED)
		message(STATUS "LTO unavailable: ${LTO_MESSAGE}")
	endif()
endif()

add_subdirectory(source)
if (BUILD_TESTS)
	add_subdirectory(tests)
endif()
if (BUILD_EXECUTABLES)
	add_subdirectory(executable)
endif()
