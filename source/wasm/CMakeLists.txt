add_executable(ausaxs_wasm wasm.cpp)
target_link_libraries(ausaxs_wasm PRIVATE ausaxs_core ausaxs_api ausaxs_math ausaxs_em ausaxs_rigidbody)
target_link_options(ausaxs_wasm PRIVATE
	--no-entry
	-fexceptions
	-sMODULARIZE=1
	-sEXPORT_NAME=createAUSAXS
	-sALLOW_MEMORY_GROWTH=1
	-sALLOW_TABLE_GROWTH=1
	-sFORCE_FILESYSTEM=1
	"-sEXPORTED_FUNCTIONS=['_malloc','_free']"
	"-sEXPORTED_RUNTIME_METHODS=['FS','ccall','cwrap','UTF8ToString','addFunction','removeFunction','HEAPF64','HEAP32','HEAPU32']"
)
set_target_properties(ausaxs_wasm PROPERTIES OUTPUT_NAME ausaxs SUFFIX ".js")
