# The backend needs nothing from the library but the ABI header, so it can also be configured on its own with
# `cmake -S source/gpu`. This spares a platform the library's own dependencies (e.g. CURL) when only the backend is wanted.
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	cmake_minimum_required(VERSION 3.15)
	project(ausaxs_gpu_backend CXX)
	set(CMAKE_CXX_STANDARD 20)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
endif()

include(GNUInstallDirs)

if (NOT ACPP_TARGETS AND NOT DEFINED ENV{ACPP_TARGETS})
	set(ACPP_TARGETS "generic" CACHE STRING "Device targets to compile the SYCL kernels for")
endif()

find_package(AdaptiveCpp CONFIG QUIET)
if (AdaptiveCpp_FOUND)
	add_library(ausaxs_gpu_sycl SHARED EXCLUDE_FROM_ALL "SYCL/SYCLBackend.cpp")
	set_target_properties(ausaxs_gpu_sycl PROPERTIES
		COMPILE_OPTIONS ""
		COMPILE_DEFINITIONS ""
		INCLUDE_DIRECTORIES ""
		LINK_LIBRARIES ""
	)
	if (WIN32)
		set_target_properties(ausaxs_gpu_sycl PROPERTIES PREFIX "") # GPULoader looks for ausaxs_gpu_sycl.dll, which a MinGW toolchain would prefix with "lib"
	endif()
	target_include_directories(ausaxs_gpu_sycl PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../include/gpu")
	target_compile_definitions(ausaxs_gpu_sycl PRIVATE AUSAXS_GPU_BUILDING_BACKEND)
	target_compile_options(ausaxs_gpu_sycl PRIVATE -O3) # the kernel is the point; never build it unoptimised
	add_sycl_to_target(TARGET ausaxs_gpu_sycl SOURCES "SYCL/SYCLBackend.cpp")
	install(TARGETS ausaxs_gpu_sycl LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} OPTIONAL)
else()
	set(_no_acpp "${CMAKE_CURRENT_BINARY_DIR}/no_adaptivecpp.cmake")
	file(WRITE ${_no_acpp}
		"message(FATAL_ERROR\n"
		"    \"the SYCL GPU backend needs an AdaptiveCpp installation, and none was found.\\n\"\n"
		"    \"Re-configure with -DAdaptiveCpp_ROOT=<AdaptiveCpp install prefix>.\"\n"
		")\n"
	)
	add_custom_target(ausaxs_gpu_sycl COMMAND ${CMAKE_COMMAND} -P ${_no_acpp} VERBATIM)
endif()
