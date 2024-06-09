#pragma once

#include <vector>
#include <string>

namespace utility {
    std::vector<std::string> split(std::string str, char delimiter);
    std::vector<std::string> split(std::string str, std::string delimiters);
    std::string join(const std::vector<std::string>& tokens, const std::string& delimiter);
    void print_warning(const std::string& text);
    void print_success(const std::string& text);
    void print_info(const std::string& text);
}