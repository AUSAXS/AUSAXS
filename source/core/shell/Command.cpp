/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <shell/Command.h>

#include <cstdio>
#include <stdexcept>
#include <string>
#include <array>

using namespace ausaxs::shell;

Command::Command() noexcept = default;
Command::Command(const std::string& cmd) : cmd(cmd) {}

void Command::set(const std::string& cmd) {
    this->cmd = cmd;
}

std::string Command::get() const {
    return cmd;
}

Command& Command::append(const shell::Argument& arg) {
    cmd += " " + arg.get();
    return *this;
}

Command& Command::append(const shell::Flag& flag) {
    cmd += " " + flag.get();
    return *this;
}

Command& Command::append(const std::string& arg) {
    cmd += " " + arg;
    return *this;
}

Command& Command::prepend(const std::string& arg) {
    cmd = arg + " " + cmd;
    return *this;
}

Command& Command::mute() {
    #ifdef __linux__
        cmd += " 2>/dev/null";
    #elif defined (_WIN32)
        cmd += " 2>NUL";
    #endif
    return *this;
}

#ifdef _WIN32
    #define popen _popen
    #define pclose _pclose
#endif
CommandResult Command::execute() const {
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