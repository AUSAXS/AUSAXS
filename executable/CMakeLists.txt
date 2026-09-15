add_executable(ausaxs "ausaxs.cpp")
target_link_libraries(ausaxs PRIVATE ausaxs_core ausaxs_api ausaxs_math ausaxs_em ausaxs_rigidbody)
set_target_properties(ausaxs PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED})