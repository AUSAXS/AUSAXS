file(GLOB_RECURSE SOURCE "${CMAKE_SOURCE_DIR}/source/md/*.cpp")
add_library(ausaxs_md OBJECT ${SOURCE})

target_include_directories(ausaxs_md PUBLIC "${CMAKE_SOURCE_DIR}/include/md")
target_link_libraries(ausaxs_md ausaxs_core)
set_target_properties(ausaxs_md PROPERTIES POSITION_INDEPENDENT_CODE 1)
