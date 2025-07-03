// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <shell/Option.h>

#include <string>

namespace ausaxs::shell {
    struct CommandResult {
        std::string out;
        int exit_code;
    };

    class Command {
        public: 
            Command() noexcept;
            Command(const std::string& cmd);

            void set(const std::string& cmd);

            std::string get() const;

            Command& append(const shell::Argument& arg);

            Command& append(const shell::Flag& flag);

            Command& append(const std::string& arg);

            Command& prepend(const std::string& arg);

            Command& mute();

            CommandResult execute() const;
            
        private: 
            std::string cmd;
    };
}