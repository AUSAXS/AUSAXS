if (NOT CMAKE_CROSSCOMPILING_EMULATOR)
	find_program(CMAKE_CROSSCOMPILING_EMULATOR node REQUIRED)
endif()

add_test(
	NAME wasm_smoke
	COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} "${CMAKE_CURRENT_SOURCE_DIR}/smoke.js" $<TARGET_FILE:ausaxs_wasm> "${CMAKE_SOURCE_DIR}/tests/files/2epe.pdb"
)

add_executable(wasm_internal internal.cpp)
target_link_libraries(wasm_internal PRIVATE ausaxs_core ausaxs_math ausaxs_em ausaxs_rigidbody ausaxs_api)
target_link_options(wasm_internal PRIVATE -fexceptions -sEXIT_RUNTIME=1)
add_test(NAME wasm_internal COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:wasm_internal>)
