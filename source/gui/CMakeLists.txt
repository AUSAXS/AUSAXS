add_library(ausaxs_gui OBJECT 
	"plotting.cpp"
	"resources.cpp"
)

target_include_directories(ausaxs_gui PUBLIC "${CMAKE_SOURCE_DIR}/include/gui")
target_link_libraries(ausaxs_gui PRIVATE ausaxs_core ausaxs_math)
set_target_properties(ausaxs_gui PROPERTIES POSITION_INDEPENDENT_CODE 1)
