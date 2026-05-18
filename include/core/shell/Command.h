// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <shell/Option.h>

#include <string>

namespace ausaxs::shell {
    /**
     * @brief The outcome of running a shell command: its captured standard output and exit code.
     */
    struct CommandResult {
        std::string out;
        int exit_code;
    };

    /**
     * @brief A shell command string that can be assembled piecewise and executed.
     */
    class Command {
        public:
            Command() noexcept;
            Command(const std::string& cmd);

            /// @brief Replace the command string with @p cmd.
            void set(const std::string& cmd);

            /// @brief Get the current command string.
            std::string get() const;

            /// @brief Append an argument (name and value) to the command.
            Command& append(const shell::Argument& arg);

            /// @brief Append a flag to the command.
            Command& append(const shell::Flag& flag);

            /// @brief Append a raw string to the command.
            Command& append(const std::string& arg);

            /// @brief Prepend a raw string to the command.
            Command& prepend(const std::string& arg);

            /// @brief Suppress the command's output by redirecting it away.
            Command& mute();

            /// @brief Run the command and return its captured output and exit code.
            CommandResult execute() const;
            
        private: 
            std::string cmd;
    };
}