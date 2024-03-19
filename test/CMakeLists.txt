set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_TEST
    "-Os -Wpedantic -Wall -Wextra"
    CACHE STRING "Flags used by the C++ compiler during test builds."
)
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

include(FetchContent)
FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG f4af9f69265d009a457aa99d1075cfba78652a66
	GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(Catch2)

# make the tests available through CTest
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)

file(GLOB_RECURSE TEST_SRC "*.cpp")

add_custom_target(tests)
set(REPORT_DIR "${CMAKE_BINARY_DIR}/test/reports")
set(TEST_BIN_DIR "${CMAKE_BINARY_DIR}/test/bin")
foreach(test ${TEST_SRC})
	cmake_path(GET test STEM testname)
	set(name "test_${testname}")
	add_executable("${name}" ${test} $<$<CONFIG:DEBUG>: ${BACKWARD_ENABLE}>)
	target_link_libraries("${name}" Catch2::Catch2WithMain)
	target_include_directories("${name}" PRIVATE "${CMAKE_SOURCE_DIR}/test")
	set_target_properties("${name}" PROPERTIES 
		RUNTIME_OUTPUT_DIRECTORY "${TEST_BIN_DIR}"
		COMPILE_FLAGS "${CMAKE_CXX_FLAGS_TEST}"
	)
	if (CMAKE_BUILD_TYPE MATCHES DEBUG)
		add_backward("${name}")
	endif()	
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
	COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_BIN_DIR}/test"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/test/files" "${TEST_BIN_DIR}/test/files"
)
