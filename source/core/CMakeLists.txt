add_library(ausaxs_core OBJECT)
set_target_properties(ausaxs_core PROPERTIES POSITION_INDEPENDENT_CODE 1)
target_include_directories(ausaxs_core PUBLIC "${CMAKE_SOURCE_DIR}/include/core")
target_link_libraries(ausaxs_core ausaxs_math CURL::libcurl)
if (DLIB) 
	target_link_libraries(ausaxs_core "dlib::dlib") 
endif()

add_library(ausaxs_static STATIC $<TARGET_OBJECTS:ausaxs_core>)
target_include_directories(ausaxs_static PUBLIC "${CMAKE_SOURCE_DIR}/include/core")
set_target_properties(ausaxs_static PROPERTIES 
	OUTPUT_NAME ausaxs_static
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)
add_library(ausaxs SHARED $<TARGET_OBJECTS:ausaxs_core>)
target_include_directories(ausaxs PUBLIC "${CMAKE_SOURCE_DIR}/include/core")
set_target_properties(ausaxs 
	PROPERTIES OUTPUT_NAME ausaxs
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)

target_link_libraries(ausaxs ausaxs_core ausaxs_math ausaxs_em ausaxs_rigidbody)
target_link_libraries(ausaxs_static ausaxs_core ausaxs_math ausaxs_em ausaxs_rigidbody)

add_subdirectory(api)
add_subdirectory(constants)
add_subdirectory(data)
add_subdirectory(dataset)
add_subdirectory(fitter)
add_subdirectory(form_factor)
add_subdirectory(grid)
add_subdirectory(hist)
add_subdirectory(hydrate)
add_subdirectory(io)
add_subdirectory(mini)
add_subdirectory(plots)
add_subdirectory(residue)
add_subdirectory(settings)
add_subdirectory(shell)
add_subdirectory(table)
add_subdirectory(utility)
