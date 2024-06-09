file(GLOB_RECURSE SOURCE "${CMAKE_SOURCE_DIR}/source/md/*.cpp")
add_library(ausaxs_gmx OBJECT ${SOURCE})

set_target_properties(ausaxs_gmx PROPERTIES
	PROPERTIES_OUTPUT_NAME ausaxs_gmx
	INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}
)
