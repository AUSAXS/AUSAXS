function(setup_compile_commands)
    option(ALLOW_UNKNOWN_ARCHITECTURE
        "Allow building on unknown CPU architectures (no -march flag applied)"
        OFF
    )

    ############################################
    ## 			Architecture mapping 		  ##
    ############################################
    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" SYS_ARCH)
    set(MARCH_FLAGS_amd64        "-march=x86-64-v3")
    set(MARCH_FLAGS_x86_64       "-march=x86-64-v3")
    set(MARCH_FLAGS_i386         "-march=i686")
    set(MARCH_FLAGS_i686         "-march=i686")
    set(MARCH_FLAGS_aarch64      "-march=armv8-a")
    set(MARCH_FLAGS_arm64        "-march=armv8-a")
    set(MARCH_FLAGS_armhf        "-march=armv7-a")
    set(MARCH_FLAGS_armv7l       "-march=armv7-a")
    set(MARCH_FLAGS_ppc64el      "-mcpu=power9")
    set(MARCH_FLAGS_ppc64le      "-mcpu=power9")
    set(MARCH_FLAGS_riscv64      "-march=rv64gc")
    set(MARCH_FLAGS_loong64      "-march=loongarch64")
    set(MARCH_FLAGS_loongarch64  "-march=loongarch64")
    set(MARCH_FLAGS_s390x        "-march=z13")

    set(MARCH_FLAG "${MARCH_FLAGS_${SYS_ARCH}}")

    if (NOT MARCH_FLAG)
        if (ALLOW_UNKNOWN_ARCHITECTURE)
            message(WARNING
                "Unknown architecture '${SYS_ARCH}'. Proceeding without -march flag.")
        else()
            message(FATAL_ERROR
                "Unknown architecture '${SYS_ARCH}'. "
                "Enable -DALLOW_UNKNOWN_ARCHITECTURE=ON to proceed anyway.")
        endif()
    endif()

    ############################################
    ## 		Common optimisation flags 		  ##
    ############################################
    set(CONSTEXPR_FLAGS
        "$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},Clang>:-fconstexpr-steps=1000000000>"
        "$<$<STREQUAL:${CMAKE_CXX_COMPILER_ID},GNU>:-fconstexpr-ops-limit=10000000000>"
    )

	set(DEBUG_FLAGS
		"$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra;-march=native>"
	)

    set(COMMON_OPT_FLAGS
        -O3
        -ffast-math
        -fno-finite-math-only
        -pipe
        ${CONSTEXPR_FLAGS}
		${DEBUG_FLAGS}
    )

    set(DEBUG_WARNINGS
        "$<$<CONFIG:DEBUG>:-g;-Wall;-Wpedantic;-Wextra>"
    )

    set(RELEASE_ARCH_FLAGS
        "$<$<CONFIG:RELEASE>:${MARCH_FLAG}>"
    )

    ############################################
    ##	   			Windows/MSVC	  		  ##
    ############################################
    if (WIN32 AND MSVC)
        add_compile_definitions(NOMINMAX;BUILD_EXPORT_DLL)
        add_compile_options(
            /fp:fast
            /constexpr:steps10000000000
            /Zm500
            /wd4267 # disable size_t --> int, unsigned int conversions
            /wd4244 # disable double --> float,int conversions
            "$<$<STREQUAL:${SYS_ARCH},x86_64>:/arch:AVX>"
            "$<$<STREQUAL:${SYS_ARCH},arm64>:/arch:armv8.0>"
        )
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded" PARENT_SCOPE)
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE PARENT_SCOPE)
        return()
    endif()

    ############################################
    ##	 			Clang/GCC 				  ##
    ############################################
    add_compile_options(
        ${COMMON_OPT_FLAGS}
        ${DEBUG_WARNINGS}
        ${RELEASE_ARCH_FLAGS}
    )

endfunction()
