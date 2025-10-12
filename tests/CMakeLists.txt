set(CMAKE_CXX_FLAGS "")
if (MSVC)
	set(CMAKE_CXX_FLAGS_TEST "/O2 /W4 /permissive- /EHsc")
	set(CMAKE_CXX_FLAGS "/EHsc")
else()
	set(CMAKE_CXX_FLAGS_TEST
		"-Os -Wpedantic -Wall -Wextra -Wno-nested-anon-types"
		CACHE STRING "Flags used by the C++ compiler during test builds."
	)
endif()

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

include(FetchContent)
FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG 914aeecfe23b1e16af6ea675a4fb5dbd5a5b8d0a
	GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(Catch2)

# make the tests available through CTest
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)

add_subdirectory(unit)
add_subdirectory(feature)

add_custom_target(tests DEPENDS unit_tests feature_tests)
