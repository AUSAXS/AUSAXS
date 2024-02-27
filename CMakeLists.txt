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
set(GUI TRUE)
set(DLIB FALSE)

include(CheckIPOSupported)
check_ipo_supported(RESULT LTO_SUPPORTED) #LTO

add_compile_definitions("$<$<CONFIG:DEBUG>:DEBUG=1;SAFE_MATH=1>")
if (WIN32)
	if (MSVC)
		add_compile_definitions("NOMINMAX")
		add_compile_options(
			/fp:fast /constexpr:steps10000000000 /Zm500 
			/wd4267 # disable size_t --> int, unsigned int conversions
			/wd4244 # disable double --> float,int conversions
		)
	elseif (MINGW)
		add_compile_options(
			-Ofast ${fprofile} -pipe -mavx
			"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},Clang>:-fconstexpr-steps=1000000000>"
			"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:-fconstexpr-ops-limit=10000000000>"
			"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra;-march=native>"
			"$<$<CONFIG:RELEASE>:-march=x86-64>"
		)
	endif()
elseif (UNIX)
#	add_compile_definitions("$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:CONSTEXPR_LOOKUP_TABLE>") # only GCC supports large constexpr arrays
	add_compile_options(
		-Ofast ${fprofile} -pipe -mavx
		"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},Clang>:-fconstexpr-steps=1000000000>"
		"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:-fconstexpr-ops-limit=10000000000>"
		"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra;-march=native>"
		"$<$<CONFIG:RELEASE>:-march=x86-64>"
	)
endif()
#add_link_options(
#	"$<$<CONFIG:DEBUG>:--coverage>"
#)

############################################
##            Dependencies                ##
############################################
include(FetchContent)
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
# when building for windows we want a single executable, so we have to include a bunch of static libraries to make curl work
if (WIN32)
	add_compile_definitions("CURL_STATICLIB")
	set(OPENSSL_USE_STATIC_LIBS TRUE)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
	set(OPENSSL_ROOT_DIR "/usr/x86_64-w64-mingw32/usr/local/lib64")
	set(CURL_LIBRARY "/usr/x86_64-w64-mingw32/usr/local/lib/libcurl.a")
	link_libraries(${CURL_LIBRARY})
	include_directories("/usr/x86_64-w64-mingw32/usr/local/include")
#	find_package(CURL REQUIRED)
	find_package(OpenSSL COMPONENTS SSL REQUIRED)
	set(LIBS -static gcc stdc++ OpenSSL::SSL winpthread ws2_32 libwldap32.a crypt32)
#	set(LIBS -static gcc stdc++ winpthread CURL::libcurl OpenSSL::SSL OpenSSL::Crypto ws2_32 libwldap32.a crypt32)

# for linux we assume the users can figure out how to install missing packages by themselves
elseif (UNIX)
	find_package(CURL REQUIRED)
	set(LIBS CURL::libcurl -static-libgcc -static-libstdc++)
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
	configure_file(${CMAKE_SOURCE_DIR}/docs/saxs.dox.in ${sim3a_Doxygen} @ONLY)
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
endif()
target_link_libraries(ausaxs ${LIBS})
target_link_libraries(ausaxs_static ${LIBS})
link_libraries(ausaxs)

### DOUBLE-COMPILED STATIC & SHARED LIBRARY APPROACH - SUPPORTS PROFILING ###
#add_library(ausaxs_static STATIC ${SOURCE_FILES})
#add_library(ausaxs SHARED ${SOURCE_FILES})
#target_include_directories(ausaxs PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
#target_include_directories(ausaxs_static PUBLIC "${CMAKE_SOURCE_DIR}/include" PRIVATE "${CMAKE_SOURCE_DIR}/source")
#target_link_libraries(ausaxs -static CLI11::CLI11 ${LIBS})
#target_link_libraries(ausaxs_static CLI11::CLI11 ${LIBS})
#link_libraries(ausaxs_static)

############################################
##             Unit Testing               ##
############################################
add_subdirectory(test)

############################################
##          Other executables             ##
############################################
add_subdirectory(executable)
