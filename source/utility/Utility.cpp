#include <utility/Utility.h>
#include <algorithm>
#include <filesystem>

#include <utility/Settings.h>
#include <utility/ConsoleColor.h>

bool utility::approx(double v1, double v2, double abs, double eps) {
    return std::abs(v1 - v2) <= std::max(abs, eps * std::max(std::abs(v1), std::abs(v2)));
}

std::string utility::remove_quotation_marks(std::string s) {
    if (s.size() > 1 && s[0] == '"' && s[s.size()-1] == '"') {
        return s.substr(1, s.size()-2);
    }
    return s;
}

std::string utility::remove_spaces(std::string s) {
    std::string::iterator end_pos = std::remove(s.begin(), s.end(), ' ');
    s.erase(end_pos, s.end());
    return s;
}

void utility::print_warning(std::string text) {
    console::print(text, console::color::red);
}

void utility::print_success(std::string text) {
    console::print(text, console::color::green);
}

void utility::print_info(std::string text) {
    console::print(text, console::color::lightblue);
}

void utility::create_directory(std::string path) {
    std::filesystem::path p(path);
    if (p.has_parent_path()) {
        std::filesystem::create_directories(p.parent_path());
    }
}

bool utility::equal(double a, double b, double c) {
    return a == b && b == c;
}

std::string utility::remove_extension(std::string path) {
    return std::filesystem::path(path).replace_extension("").string();
}

std::string utility::stem_append(std::string path, std::string s) {
    std::filesystem::path p(path);
    return p.parent_path().string() + "/" + p.stem().string() + s + p.extension().string();
}

std::string utility::stem(std::string path) {
    return std::filesystem::path(path).stem().string();    
}

std::vector<std::string> utility::split(std::string str, char delimiter) {
    std::string token;
    std::stringstream ss(str);
    std::vector<std::string> tokens;
    while(std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::string> utility::split(std::string str, std::string delimiters) {
    std::vector<std::string> tokens;

    auto is_delimiter = [&delimiters] (char c) {
        return delimiters.find(c) != std::string::npos;
    };

    // skip leading delimiters
    unsigned int start = 0;
    for (; start < str.size(); start++) {
        if (!is_delimiter(str[start])) {
            break;
        }
    }

    // iterate through the rest of the string
    for (unsigned int i = start; i < str.size(); i++) {
        if (!is_delimiter(str[i])) {
            continue;
        }

        // add token to vector
        tokens.push_back(str.substr(start, i-start));
        start = ++i;

        // skip consecutive delimiters
        for (; start < str.size(); start++) {
            if (!is_delimiter(str[start])) {
                break;
            }
        }
        i = start;
    }

    // add last token to vector
    if (start < str.size()) {
        tokens.push_back(str.substr(start));
    }
    return tokens;
}

std::string utility::join(std::vector<std::string> v, std::string separator) {
    std::string s;
    for (unsigned int i = 0; i < v.size(); i++) {
        s += v[i];
        if (i != v.size()-1) {
            s += separator;
        }
    }
    return s;
}

std::string utility::remove_all(std::string s, std::string remove) {
    std::string new_s;
    for (auto c : s) {
        if (remove.find(c) == std::string::npos) {
            new_s += c;
        }
    }
    return new_s;
}

std::string utility::to_lowercase(std::string s) {
    std::string new_s;
    for (auto c : s) {
        new_s += std::tolower(c);
    }
    return new_s;
}

std::string utility::uid() {
    static unsigned int i = 0;
    return std::to_string(i++);
}

std::string utility::uid(std::string s) {return s + uid();}

std::ostream& utility::detail::operator<<(std::ostream& os, const __dummy& obj) {
    os << obj.s;
    return os;
}

utility::detail::__dummy utility::fixedwidth(double number, unsigned int width) {
    std::string s = std::to_string(number);
    std::string o;
    for (unsigned int i = 0; i < width; i++) {
        if (i < s.size()) {
            o += s[i];
        } else {
            o += ' ';
        }
    }

    // check how lossy the conversion was
    double d = std::stod(o);
    if (!approx(d, number, 1e-3)) {
        print_warning("Fixed-width conversion of " + std::to_string(number) + " to " + o + " is lossy.");
    }
    
    return {o};
} 