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

add_subdirectory(unit)
add_subdirectory(feature)
add_subdirectory(benchmark)

add_custom_target(tests DEPENDS unit_tests feature_tests)
