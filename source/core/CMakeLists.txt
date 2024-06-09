set(mainlib "objlib")
add_library(${mainlib} OBJECT)
set_target_properties(${mainlib} PROPERTIES POSITION_INDEPENDENT_CODE 1)

add_library(ausaxs_static STATIC $<TARGET_OBJECTS:objlib>)
add_library(ausaxs SHARED $<TARGET_OBJECTS:objlib>)
set_target_properties(ausaxs_static PROPERTIES 
	OUTPUT_NAME ausaxs_static
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)
set_target_properties(ausaxs 
	PROPERTIES OUTPUT_NAME ausaxs
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)

target_include_directories(${mainlib} PUBLIC "${CMAKE_SOURCE_DIR}/include/core")
if (DLIB)
	target_link_libraries(${mainlib} dlib::dlib)
endif()
target_link_libraries(ausaxs ${LIBS})
target_link_libraries(ausaxs_static ${LIBS})

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
