// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <string>

#ifdef _WIN32
    #include <process.h>
    #define AUSAXS_TEST_GETPID _getpid
#else
    #include <unistd.h>
    #define AUSAXS_TEST_GETPID getpid
#endif

namespace ausaxs::test {
    /**
     * @brief A scratch file path unique to this process.
     *
     * catch_discover_tests registers every TEST_CASE as its own CTest test, and CI runs them with `ctest --parallel`. Several processes of the *same* test
     * binary are therefore alive at once, each with its own copy of any function-local counter. A path built from `prefix + counter` alone collides across
     * them: one process truncates and rewrites the file while another is parsing it, so a test silently runs someone else's input. Mixing the pid in makes
     * the path unique per process, which is what the counter was assumed to provide.
     *
     * @param prefix File name prefix, without a directory.
     * @param counter A per-process counter distinguishing the paths within one process.
     * @param extension File extension, including the dot.
     */
    inline std::string temp_path(const std::string& prefix, int counter, const std::string& extension) {
        return "/tmp/" + prefix + "_" + std::to_string(AUSAXS_TEST_GETPID()) + "_" + std::to_string(counter) + extension;
    }
}

#undef AUSAXS_TEST_GETPID
