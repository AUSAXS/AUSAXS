add_library(ausaxs_api OBJECT 
    "ObjectStorage.cpp"

    "api_helper.cpp"
    "api_sasview.cpp"

    "pyausaxs/api_data.cpp"
    "pyausaxs/api_fit.cpp"
    "pyausaxs/api_iterative_fit.cpp"
    "pyausaxs/api_molecule.cpp"
    "pyausaxs/api_pdb.cpp"
    "pyausaxs/api_settings.cpp"
)

target_include_directories(ausaxs_api PUBLIC "${CMAKE_SOURCE_DIR}/include/api")
target_link_libraries(ausaxs_api ausaxs_core ausaxs_rigidbody)
set_target_properties(ausaxs_api PROPERTIES POSITION_INDEPENDENT_CODE 1)
