file(GLOB_RECURSE TEST_SRC "*.cpp")

add_custom_target(feature_tests)
set(REPORT_DIR "${CMAKE_BINARY_DIR}/tests/reports")
set(TEST_BIN_DIR "${CMAKE_BINARY_DIR}/tests/feature/bin")
foreach(test ${TEST_SRC})
	cmake_path(GET test STEM testname)
	set(name "ftest_${testname}")
	add_executable("${name}" ${test})
	target_link_libraries("${name}" Catch2::Catch2WithMain ausaxs "$<$<CONFIG:Debug>:Backward::Interface>")
	target_include_directories("${name}" PRIVATE "${CMAKE_SOURCE_DIR}/tests/feature")
	set_target_properties("${name}" PROPERTIES 
		RUNTIME_OUTPUT_DIRECTORY "${TEST_BIN_DIR}"
		COMPILE_FLAGS "${CMAKE_CXX_FLAGS_TEST}"
	)
	add_dependencies(feature_tests "${name}")
	catch_discover_tests("${name}" 
		TEST_SPEC "~[slow] ~[broken] ~[manual]" 
		WORKING_DIRECTORY "${TEST_BIN_DIR}"
		REPORTER junit
		OUTPUT_DIR ${REPORT_DIR}
		OUTPUT_SUFFIX ".xml"
	)
endforeach()

# since the tests are run in the TEST_BIN_DIR directory, we have to create a symbolic link from there to the test files
add_custom_command(TARGET feature_tests POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory "${REPORT_DIR}"
	COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_BIN_DIR}/tests"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/tests/files" "${TEST_BIN_DIR}/tests/files"
)
