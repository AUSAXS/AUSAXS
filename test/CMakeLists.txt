set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE TEST)

include(FetchContent)
FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG f4af9f69265d009a457aa99d1075cfba78652a66
	GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(Catch2)

file(GLOB_RECURSE 
	TEST_SRC
	"*.cpp"
)

# make the tests available through cmake
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)

add_custom_target(tests)
foreach(test ${TEST_SRC})
	cmake_path(GET test STEM testname)
	set(name "test_${testname}")
	add_executable("${name}" ${test} ${ENABLE_DEBUG})
	target_link_libraries("${name}" Catch2::Catch2WithMain)
	set_target_properties("${name}" PROPERTIES 
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test/bin"
	)
	if (CMAKE_BUILD_TYPE MATCHES DEBUG)
		add_backward("${name}")
	endif()	
	add_dependencies(tests "${name}")
	catch_discover_tests("${name}" TEST_SPEC "~[slow] ~[broken] ~[manual]" WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/..")
endforeach()
