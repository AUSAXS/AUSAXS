#pragma once

#include <md/shell/Option.h>

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

namespace ausaxs::shell {
    struct CommandResult {
        std::string out;
        int exit_code;
    };

    class Command {
        public: 
            Command() = default;
            Command(const std::string& cmd) : cmd(cmd) {}

            void set(const std::string& cmd) {
                this->cmd = cmd;
            }

            std::string get() const {
                return cmd;
            }

            Command& append(shell::Argument arg) {
                cmd += " " + arg.get();
                return *this;
            }

            Command& append(shell::Flag flag) {
                cmd += " " + flag.get();
                return *this;
            }

            Command& append(std::string arg) {
                cmd += " " + arg;
                return *this;
            }

            Command& prepend(std::string arg) {
                cmd = arg + " " + cmd;
                return *this;
            }

            Command& mute() {
                #ifdef __linux__
                    cmd += " 2>/dev/null";
                #elif _WIN32
                    cmd += " 2>NUL";
                #endif
                return *this;
            }

            CommandResult execute() const {
                std::array<char, 128> buffer;
                std::string result;
                FILE* pipe = popen(cmd.data(), "r");
                if (pipe == nullptr) {
                    throw std::runtime_error("popen() failed!");
                }
                while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
                    result += buffer.data();
                }
                int exit_code = pclose(pipe);
                return {std::move(result), exit_code};
            }
            
        private: 
            std::string cmd;
    };
}