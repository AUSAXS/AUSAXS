cmake_minimum_required(VERSION 3.14)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE DEBUG)

project(SAXS VERSION 0.0)
#set(fprofile "-fprofile-generate=${CMAKE_BINARY_DIR}/profiling")
#set(fprofile "-fprofile-use=${CMAKE_BINARY_DIR}/profiling" "-fprofile-correction")
#add_link_options(-lgcov --coverage)
set(fprofile "")

include(CheckIPOSupported)
check_ipo_supported(RESULT LTO_SUPPORTED) #LTO

add_compile_definitions("$<$<CONFIG:DEBUG>:DEBUG=1;SAFE_MATH=1>")#;_GLIBCXX_DEBUG>")
if (MSVC)
	add_compile_options(
		/O2 /ffinite-math /constexpr:steps10000000000
	)
elseif (UNIX)
	add_compile_options(
		-Ofast ${fprofile} -pipe -mavx -fconstexpr-ops-limit=10000000000
		"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra;-march=native>"
		"$<$<CONFIG:RELEASE>:-march=native>"#-march=x86-64>"
	)
endif()
#add_link_options(
#	"$<$<CONFIG:DEBUG>:--coverage>"
#)

############################################
##            Dependencies                ##
############################################
include(FetchContent)
#set(FETCHCONTENT_UPDATES_DISCONNECT ON)
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

FetchContent_Declare(
	CLI11
	GIT_REPOSITORY https://github.com/CLIUtils/CLI11
	GIT_TAG 792d89286788acac125e0487f8dbde88035f7422
	GIT_PROGRESS TRUE
)

FetchContent_Declare(
	backward
	GIT_REPOSITORY https://github.com/bombela/backward-cpp
	GIT_TAG v1.6
)

FetchContent_Declare(
	BSThreadPool
	GIT_REPOSITORY https://github.com/klytje/thread-pool
)

FetchContent_Declare(
	gcem
	GIT_REPOSITORY https://github.com/kthohr/gcem.git
)

set(CMAKE_CXX_STANDARD 20) # continue with C++20
FetchContent_MakeAvailable(CLI11 BSThreadPool GCEM)
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    FetchContent_MakeAvailable(backward)
endif()
include_directories(${CMAKE_BINARY_DIR}/_deps/bsthreadpool-src ${gcem_SOURCE_DIR}/include)
#add_compile_options("$<$<CONFIG:DEBUG>:-fopt-info-loop;-fopt-info-vec>")

############################################
##           Find and link CURL           ##
############################################
find_package(CURL REQUIRED)
find_package(OpenSSL COMPONENTS SSL REQUIRED)

# when building for windows we want a single executable, so we have to include a bunch of static libraries to make curl work
if (WIN32)
	set(LIBS -static gcc stdc++ winpthread CURL::libcurl OpenSSL::SSL OpenSSL::Crypto ws2_32 libwldap32.a crypt32)

# for linux we assume the users can figure out how to install missing packages by themselves
elseif (UNIX)
	set(LIBS CURL::libcurl OpenSSL::SSL)
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

### SINGLE-COMPILED OBJECT LIBRARY APPROACH ###
#add_library(objlib OBJECT ${SOURCE_FILES})
#set_target_properties(objlib PROPERTIES POSITION_INDEPENDENT_CODE 1)

#add_library(saxs_static STATIC $<TARGET_OBJECTS:objlib>)
#add_library(saxs SHARED $<TARGET_OBJECTS:objlib>)
#set_target_properties(saxs_static PROPERTIES 
#	OUTPUT_NAME saxs_static
#	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
#)
#set_target_properties(saxs 
#	PROPERTIES OUTPUT_NAME saxs
#	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
#)

#target_include_directories(objlib PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
#target_include_directories(saxs PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
#target_include_directories(saxs_static PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
#target_link_libraries(objlib CLI11::CLI11 dlib::dlib ${LIBS})
#target_link_libraries(saxs CLI11::CLI11 dlib::dlib ${LIBS})
#target_link_libraries(saxs_static CLI11::CLI11 dlib::dlib ${LIBS})
#link_libraries(saxs_static)

### DOUBLE-COMPILED STATIC & SHARED LIBRARY APPROACH - SUPPORTS PROFILING ###
add_library(saxs_static STATIC ${SOURCE_FILES})
add_library(saxs SHARED ${SOURCE_FILES})
target_include_directories(saxs PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
target_include_directories(saxs_static PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
target_link_libraries(saxs CLI11::CLI11 dlib::dlib ${LIBS})
target_link_libraries(saxs_static CLI11::CLI11 dlib::dlib ${LIBS})
link_libraries(saxs_static)

############################################
##             Unit Testing               ##
############################################
add_subdirectory(test)

############################################
##          Other executables             ##
############################################
add_subdirectory(executable)
