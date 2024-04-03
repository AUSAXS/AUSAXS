set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_minimum_required(VERSION 3.15)
project(AUSAXS VERSION 1.0)
option(GUI "Enable GUI executables" ON)
option(DLIB "Download and use the dlib minimizers" ON)
option(BUILD_PLOT_EXE "Compile the plotting utility as an executable for the current platform" FALSE)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_compile_definitions("$<$<CONFIG:DEBUG>:DEBUG=1;SAFE_MATH=1>")
if (WIN32)
	if (MSVC)
		add_compile_definitions("NOMINMAX")
		add_compile_options(
			/fp:fast /constexpr:steps10000000000 /Zm500 
			/wd4267 # disable size_t --> int, unsigned int conversions
			/wd4244 # disable double --> float,int conversions
		)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
	elseif (MINGW)
		add_compile_options(
			-Ofast -pipe -mavx
			"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},Clang>:-fconstexpr-steps=1000000000>"
			"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:-fconstexpr-ops-limit=10000000000>"
			"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra;-march=native>"
			"$<$<CONFIG:RELEASE>:-march=x86-64>"
		)
	endif()
elseif (UNIX)
	add_compile_options(
		-Ofast -pipe -mavx
		"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},Clang>:-fconstexpr-steps=1000000000>"
		"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:-fconstexpr-ops-limit=10000000000>"
		"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra;-march=native>"
		"$<$<CONFIG:RELEASE>:-march=x86-64>"
	)
endif()

############################################
##            Dependencies                ##
############################################
include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
if (DLIB)
	FetchContent_Declare(
		dlib
		GIT_REPOSITORY https://github.com/davisking/dlib
		GIT_TAG v19.24.2
		GIT_PROGRESS TRUE
	)

	# tell dlib not to link with various unnecessary libraries
	SET(DLIB_NO_GUI_SUPPORT TRUE)
	SET(DLIB_JPEG_SUPPORT OFF)
	SET(DLIB_LINK_WITH_SQLITE3 OFF)
	SET(DLIB_USE_BLAS OFF)
	SET(DLIB_USE_LAPACK OFF)
	SET(DLIB_USE_CUDA OFF)
	SET(DLIB_PNG_SUPPORT OFF)
	SET(DLIB_GIF_SUPPORT OFF)
	SET(DLIB_WEBP_SUPPORT OFF)
	SET(DLIB_USE_FFTW OFF)
	SET(DLIB_USE_MKL_FFT OFF)
	SET(DLIB_USE_FFMPEG OFF)
	set(CMAKE_CXX_STANDARD 17) # dlib must be compiled with C++17
	FetchContent_MakeAvailable(dlib)
	set(CMAKE_CXX_STANDARD 20) # continue with C++20	
	add_compile_definitions("DLIB_AVAILABLE")
endif()

FetchContent_Declare(
	backward
	GIT_REPOSITORY https://github.com/bombela/backward-cpp
	GIT_TAG v1.6
)

FetchContent_Declare(
	thread_pool
	GIT_REPOSITORY https://github.com/bshoshany/thread-pool
)

FetchContent_Declare(
	gcem
	GIT_REPOSITORY https://github.com/kthohr/gcem
)

FetchContent_MakeAvailable(thread_pool GCEM backward)
include_directories(${thread_pool_SOURCE_DIR}/include ${gcem_SOURCE_DIR}/include ${backward_SOURCE_DIR})

############################################
##           Find and link CURL           ##
############################################
if (WIN32)
	add_compile_definitions("CURL_STATICLIB")
	find_package(CURL REQUIRED)
	find_package(OpenSSL COMPONENTS SSL REQUIRED)
	set(LIBS -static CURL::libcurl)
	include_directories(${CURL_INCLUDE_DIRS}) # find_package apparently does not always include the headers
elseif (UNIX)
	find_package(CURL REQUIRED)
	set(LIBS CURL::libcurl -static-libgcc -static-libstdc++)
	if (CMAKE_BUILD_TYPE MATCHES DEBUG)
		link_libraries(-lbfd -ldl)
	endif()
endif()


############################################
##                Doxygen                 ##
############################################
find_package(Doxygen)
if(DOXYGEN_FOUND)
	set(sim3a_Doxygen "${CMAKE_BINARY_DIR}/saxs.dox")
	configure_file(${CMAKE_SOURCE_DIR}/scripts/ausaxs.dox.in ${sim3a_Doxygen} @ONLY)
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
include(CheckIPOSupported)
check_ipo_supported(RESULT LTO_SUPPORTED)

file(GLOB_RECURSE HEADER_FILES "${CMAKE_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE SOURCE_FILES "${CMAKE_SOURCE_DIR}/source/*.cpp")

add_library(objlib OBJECT ${SOURCE_FILES})
set_target_properties(objlib PROPERTIES POSITION_INDEPENDENT_CODE 1)

add_library(ausaxs_static STATIC $<TARGET_OBJECTS:objlib>)
add_library(ausaxs SHARED $<TARGET_OBJECTS:objlib>)
set_target_properties(ausaxs_static PROPERTIES 
	OUTPUT_NAME ausaxs_static
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)
set_target_properties(ausaxs 
	PROPERTIES OUTPUT_NAME ausaxs
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)

target_include_directories(objlib PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
target_include_directories(ausaxs PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
target_include_directories(ausaxs_static PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
if (DLIB)
	target_link_libraries(objlib dlib::dlib)
	target_link_libraries(ausaxs dlib::dlib)
	target_link_libraries(ausaxs_static dlib::dlib)
endif()
target_link_libraries(ausaxs ${LIBS})
target_link_libraries(ausaxs_static ${LIBS})

add_subdirectory(test)
add_subdirectory(executable)
