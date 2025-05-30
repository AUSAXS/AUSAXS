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

file(GLOB_RECURSE TEST_SRC "*.cpp")

add_custom_target(tests)
set(REPORT_DIR "${CMAKE_BINARY_DIR}/tests/reports")
set(TEST_BIN_DIR "${CMAKE_BINARY_DIR}/tests/bin")
foreach(test ${TEST_SRC})
	cmake_path(GET test STEM testname)
	set(name "test_${testname}")
	add_executable("${name}" ${test})
	target_link_libraries("${name}" Catch2::Catch2WithMain ausaxs "$<$<CONFIG:Debug>:Backward::Interface>")
	target_include_directories("${name}" PRIVATE "${CMAKE_SOURCE_DIR}/tests")
	set_target_properties("${name}" PROPERTIES 
		RUNTIME_OUTPUT_DIRECTORY "${TEST_BIN_DIR}"
		COMPILE_FLAGS "${CMAKE_CXX_FLAGS_TEST}"
	)
	add_dependencies(tests "${name}")
	catch_discover_tests("${name}" 
		TEST_SPEC "~[slow] ~[broken] ~[manual]" 
		WORKING_DIRECTORY "${TEST_BIN_DIR}"
		REPORTER junit
		OUTPUT_DIR ${REPORT_DIR}
		OUTPUT_SUFFIX ".xml"
	)
endforeach()

# since the tests are run in the TEST_BIN_DIR directory, we have to create a symbolic link from there to the test files
add_custom_command(TARGET tests POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory "${REPORT_DIR}"
	COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_BIN_DIR}/tests"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/tests/files" "${TEST_BIN_DIR}/tests/files"
)
