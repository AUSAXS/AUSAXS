set(CMAKE_CXX_FLAGS "")
if (MSVC)
	set(CMAKE_CXX_FLAGS "/EHsc")
	set(AUSAXS_TEST_FLAGS /O2 /W4 /permissive-)
else()
	set(AUSAXS_TEST_FLAGS -Os -Wpedantic -Wall -Wextra -Wno-nested-anon-types)
endif()

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

if(USE_SYSTEM_CATCH)
	find_package(Catch2 REQUIRED)
else()
	include(FetchContent)
	FetchContent_Declare(
		Catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG 914aeecfe23b1e16af6ea675a4fb5dbd5a5b8d0a
		GIT_PROGRESS TRUE
	)
	FetchContent_MakeAvailable(Catch2)
endif()

# make the tests available through CTest
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)

# linked into every unit and feature test; see the file for what it does
add_library(ausaxs_test_support OBJECT support/gpu_listener.cpp)
target_link_libraries(ausaxs_test_support PUBLIC Catch2::Catch2 libausaxs)
target_compile_options(ausaxs_test_support PRIVATE ${AUSAXS_TEST_FLAGS})

add_subdirectory(unit)
add_subdirectory(feature)
add_subdirectory(benchmark)

add_custom_target(tests DEPENDS unit_tests feature_tests)
