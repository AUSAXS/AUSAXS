#include <utility/Utility.h>
#include <algorithm>
#include <filesystem>

#include <settings.h>

bool utility::approx(double v1, double v2, double abs, double eps) {
    if (v1-abs > v2*(1+eps)) {return false;}
    if (v1+abs < v2*(1-eps)) {return false;}
    return true;
}

// std::string remove_spaces(std::string s) {
//     std::string::iterator end_pos = std::remove(s.begin(), s.end(), ' ');
//     s.erase(end_pos, s.end());
//     return s;
// }

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