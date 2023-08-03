file(GLOB_RECURSE EXEC_SOURCE "*.cpp")

foreach(exec ${EXEC_SOURCE})
	cmake_path(GET exec STEM name)
	add_executable("${name}" ${exec} $<$<CONFIG:DEBUG>: ${BACKWARD_ENABLE}>)
	set_target_properties("${name}" PROPERTIES 
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
	)
	if (CMAKE_BUILD_TYPE MATCHES DEBUG)
		add_backward("${name}")
	endif()	
endforeach()
