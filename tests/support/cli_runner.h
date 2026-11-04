// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <support/temp_file.h>

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>

#ifndef AUSAXS_EXECUTABLE
    #error "cli_runner.h needs AUSAXS_EXECUTABLE, which tests/feature/CMakeLists.txt defines for tests under tests/feature/cli/ only."
#endif

#ifndef _WIN32
    #include <sys/wait.h>
#endif

namespace ausaxs::test {
    /**
     * @brief What a user sees from one invocation of the ausaxs executable.
     */
    struct CLIResult {
        // the process exit code, with a signal reported the way a shell reports it (128 + signal), so an abort is 134
        int exit_code = 0;

        // standard output and standard error, interleaved in the order they were written
        std::string output;

        // console:: wraps every message in ANSI colour codes, so a whole printed line never matches literally - the message text between them does.
        // Match on a distinctive fragment, not on the line.
        [[nodiscard]] bool contains(std::string_view text) const {return output.find(text) != std::string::npos;}
    };

    /**
     * @brief Run the ausaxs executable with the given arguments and capture its exit code and terminal output.
     *
     * This is the only way to reach the CLI layer itself: argument handling, the exit code and what reaches the terminal are properties of the
     * executable, and a test calling cli_saxs_fitter() directly exercises none of them. Both streams are redirected into one file through the shell,
     * which is what makes the output ordering match what a user would see.
     *
     * @param arguments The arguments to pass, as one string, e.g. "fit tests/files/2epe.pdb tests/files/2epe.dat -o temp/out/".
     */
    inline CLIResult run_cli(std::string_view arguments) {
        TempFile log(".log");
        int status = std::system((std::string("\"") + AUSAXS_EXECUTABLE + "\" " + std::string(arguments) + " > \"" + log.path() + "\" 2>&1").c_str());

        CLIResult result;
#ifdef _WIN32
        result.exit_code = status;
#else
        // std::system hands back a wait status, not an exit code. Decode it as a shell would, so a crash shows up as 128 + signal (an abort is 134)
        // rather than as some arbitrary large number.
        if (WIFEXITED(status)) {result.exit_code = WEXITSTATUS(status);}
        else if (WIFSIGNALED(status)) {result.exit_code = 128 + WTERMSIG(status);}
        else {result.exit_code = status;}
#endif

        std::ifstream input(log.path(), std::ios::binary);
        result.output.assign(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
        return result;
    }
}
