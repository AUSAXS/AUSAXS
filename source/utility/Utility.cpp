#include <utility/Utility.h>
#include <algorithm>
#include <filesystem>

#include <utility/Settings.h>

bool utility::approx(double v1, double v2, double abs, double eps) {
    if (v1-abs > v2*(1+eps)) {return false;}
    if (v1+abs < v2*(1-eps)) {return false;}
    return true;
}

std::string utility::remove_spaces(std::string s) {
    std::string::iterator end_pos = std::remove(s.begin(), s.end(), ' ');
    s.erase(end_pos, s.end());
    return s;
}

void utility::print_warning(std::string text) {
    std::cout << "\033[1;31m" << text << "\033[0m" << std::endl;
}

void utility::create_directories(std::string& path) {
    std::filesystem::path p(path);
    if (p.has_parent_path()) {
        std::filesystem::create_directories(p.parent_path());
    }

    if (!p.has_extension()) {
        path += "." + setting::figures::format;
    }
}

std::string utility::remove_extension(std::string path) {
    return std::filesystem::path(path).replace_extension("");
}

std::string utility::stem_append(std::string path, std::string s) {
    std::filesystem::path p(path);
    return p.parent_path().string() + "/" + p.stem().string() + s + p.extension().string();
}

template<>
std::string utility::extract_number<std::string>(std::string s) {
    unsigned int start = 0;
    while (!std::isdigit(s[start]) && start != s.size()) {start++;}
    unsigned int end = start;
    while ((std::isdigit(s[end]) || s[end] == '.') && end != s.size()) {end++;}
    while (end > 0 && s[end-1] == '.') {end--;}
    return s.substr(start, end-start);
}

std::string utility::uid() {
    static unsigned int i = 0;
    return std::to_string(i++);
}

std::string utility::uid(std::string s) {return s + uid();}