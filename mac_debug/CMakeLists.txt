set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_minimum_required(VERSION 3.15)
project(AUSAXS VERSION 1.0)
option(GUI "Enable GUI executables" ON)
option(DLIB "Download and use the dlib minimizers" ON)
option(BUILD_PLOT_EXE "Compile the plotting utility as an executable for the current platform" FALSE)
set(ARCH "native" CACHE STRING "Target architecture. Default: native")

add_compile_definitions("$<$<CONFIG:DEBUG>:DEBUG=1;SAFE_MATH=1>")
add_compile_options(
	-Ofast -pipe
	"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},Clang>:-fconstexpr-steps=1000000000>"
	"$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:-fconstexpr-ops-limit=10000000000>"
	"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra;-march=native>"
	"$<$<AND:$<CONFIG:RELEASE>,$<STREQUAL:${ARCH},native>>:-march=native>"
	"$<$<AND:$<CONFIG:RELEASE>,$<STREQUAL:${ARCH},x86-64>>:-march=x86-64-v3 -mavx>"
	"$<$<AND:$<CONFIG:RELEASE>,$<STREQUAL:${ARCH},arm64>>:-march=armv8-a>"
)

find_package(CURL REQUIRED)
set(LIBS -static CURL::libcurl)

add_executable(hello_world main.cpp)
target_link_libraries(hello_world ${LIBS})
