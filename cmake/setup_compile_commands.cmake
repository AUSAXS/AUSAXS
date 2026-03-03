function(setup_compile_commands)
    option(ALLOW_UNKNOWN_ARCHITECTURE
        "Allow building on unknown CPU architectures (no -march flag applied)"
        OFF
    )

    ############################################
    ##          Architecture mapping          ##
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

    if (NOT DEFINED MARCH_FLAGS_${SYS_ARCH})
        if (ALLOW_UNKNOWN_ARCHITECTURE)
            message(WARNING
                "Unknown architecture '${SYS_ARCH}'. Proceeding without -march flag.")
        else()
            message(FATAL_ERROR
                "Unknown architecture '${SYS_ARCH}'. "
                "Enable -DALLOW_UNKNOWN_ARCHITECTURE=ON to proceed anyway.")
        endif()
    endif()
    set(MARCH_FLAG "${MARCH_FLAGS_${SYS_ARCH}}")

    ############################################
    ##    Per-compiler optimisation flags     ##
    ############################################
    set(CompilerFlags "")

    if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        # MSVC uses /arch: flags below; -march is not applicable
        add_compile_definitions(NOMINMAX BUILD_EXPORT_DLL)
        list(APPEND CompilerFlags
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

    else()
        list(APPEND CompilerFlags
            -O3
            -ffast-math
            -fno-finite-math-only
            -pipe
            "$<$<CONFIG:Debug>:-g>"
            "$<$<CONFIG:Debug>:-Wall>"
            "$<$<CONFIG:Debug>:-Wpedantic>"
            "$<$<CONFIG:Debug>:-Wextra>"
        )

        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            list(APPEND CompilerFlags "-fconstexpr-ops-limit=10000000000")

        elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            list(APPEND CompilerFlags "-fconstexpr-steps=1000000000")

        elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
            list(APPEND CompilerFlags "-fconstexpr-steps=1000000000")
            # -march is not supported on arm64 Apple targets; the SDK/toolchain handles targeting.
            # Use CMAKE_OSX_ARCHITECTURES when explicitly set (CI), otherwise fall back to the
            # host processor. Suppress -march whenever the target includes arm64 (covers both
            # single arm64 and universal arm64;x86_64 builds).
            if (CMAKE_OSX_ARCHITECTURES)
                set(_MARCH_TARGET "${CMAKE_OSX_ARCHITECTURES}")
            else()
                set(_MARCH_TARGET "${SYS_ARCH}")
            endif()
            if (_MARCH_TARGET MATCHES "arm64")
                set(MARCH_FLAG "")
            endif()

        endif()

        list(APPEND CompilerFlags ${MARCH_FLAG})
    endif()

    add_compile_options(${CompilerFlags})

endfunction()
