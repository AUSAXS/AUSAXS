#include <utility/Settings.h>
#include <utility/Utility.h>

#include <fstream>
#include <filesystem>

bool is_comment_char(char c) {
    switch (c) {
        case '/': [[fallthrough]];
        case '#': [[fallthrough]];
        case '%': return true;
        default: return false;
    }
}

void set(std::string opt, std::vector<std::string> val) {
    for (const auto& e : setting::detail::options) {
        for (const auto& alias : e->aliases) {
            if (alias == opt) {
                e->set(val);
                return;
            }
        }
    }

    throw except::unexpected("Settings::set: Option \"" + opt + "\" not found.");
}

void setting::read(std::string path) {
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("Settings::read: Could not open setup file.");}

    auto split_tokens = [] (std::string line) {
        unsigned int start = 0;

        // find end of first word
        while (start < line.size() && line[start] != ' ') {
            start++;
            continue;
        }
        if (start == line.size()) {return std::pair<std::string, std::vector<std::string>>("empty", {});} // line consists only of spaces
        std::string first = line.substr(0, start);

        unsigned int end = start;
        std::vector<std::string> vals;
        while (end != line.size()) {
            // find start of word, skipping any amount of spacing
            while (start < line.size() && line[start] == ' ') {
                start++;
                continue;
            }
            if (is_comment_char(line[start])) {break;} // stop if we find a comment character
            end = start;

            // find end of word
            while (line[end] != ' ' && end < line.size()) {
                end++;
                continue;
            }
            vals.push_back(line.substr(start, end-start));
            start = end;
        }
        if (vals.empty()) {throw except::parse_error("Settings::read: Could not find matching value for option \"" + first + "\"");}
        return std::pair(first, vals);
    };

    std::string line; 
    while (getline(input, line)) {
        if (is_comment_char(line[0])) {continue;} // skip comments
        if (line.empty()) {continue;}             // skip empty lines
        auto[first, second] = split_tokens(line);
        if (first == "empty") {continue;}         // skip empty lines
        set(first, second);
    }
}

void setting::write(std::string path) {
    utility::create_directory(path);
    std::ofstream output(path);
    if (!output.is_open()) {throw std::ios_base::failure("Settings::read: Could not open setup file.");}

    // write settings
    output << "### Auto-generated settings file ###\n";
    for (const auto& e : detail::options) {
        output << e->aliases[0] << " " << e->get() << "\n";
    }
    output.close();
}

bool setting::discover(std::string path) {
    if (path[path.size()-1] != '/') {path += "/";}
    std::vector<std::string> valid_names = {"settings", "setting", "setup", "config"};
    for (const auto& e : valid_names) {
        if (std::filesystem::exists(path + e + ".txt")) {
            std::cout << "Using discovered settings file at " << path + e + ".txt" << std::endl;
            setting::read(path + e + ".txt");
            return true;
        }
    }
    return false;
}

template<>
std::string setting::detail::SmartOption<std::vector<std::string>>::get() const {
    std::string str;
    std::for_each(setting.begin(), setting.end(), [&str] (std::string s) {str += s + " ";});
    return str;
}

template<>
std::string setting::detail::SmartOption<std::vector<double>>::get() const {
    std::string str;
    std::for_each(setting.begin(), setting.end(), [&str] (double s) {str += std::to_string(s) + " ";});
    return str;
}

template<>
std::string setting::detail::SmartOption<std::string>::get() const {return setting;}

template<>
std::string setting::detail::SmartOption<double>::get() const {return std::to_string(setting);}

template<>
std::string setting::detail::SmartOption<int>::get() const {return std::to_string(setting);}

template<>
std::string setting::detail::SmartOption<unsigned int>::get() const {return std::to_string(setting);}

template<>
std::string setting::detail::SmartOption<bool>::get() const {return std::to_string(setting);}

template<>
void setting::detail::SmartOption<std::string>::set(std::vector<std::string> str) const {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option received too many values.");}
    setting = str[0];
}

template<>
void setting::detail::SmartOption<bool>::set(std::vector<std::string> str) const {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option received too many values.");}

    if (str[0] == "true" || str[0] == "TRUE" || str[0] == "1") {setting = true; return;}
    else if (str[0] == "false" || str[0] == "FALSE" || str[0] == "0") {setting = false; return;}
    throw except::parse_error("Settings::parse_bool: Expected boolean std::string, but got \"" + str[0] + "\".");
}

template<>
void setting::detail::SmartOption<double>::set(std::vector<std::string> str) const {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option received too many values.");}
    setting = std::stod(str[0]);
}

template<>
void setting::detail::SmartOption<int>::set(std::vector<std::string> str) const {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option received too many values.");}
    setting = std::stoi(str[0]); 
}

template<>
void setting::detail::SmartOption<unsigned int>::set(std::vector<std::string> str) const {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option received too many values.");}
    setting = std::stoi(str[0]);
}

template<>
void setting::detail::SmartOption<std::vector<std::string>>::set(std::vector<std::string> str) const {
    setting = str;
}

template<>
void setting::detail::SmartOption<std::vector<double>>::set(std::vector<std::string> str) const {
    std::vector<double> new_val;
    for (auto& s : str) {
        if (s.empty() || s == " ") {continue;}
        new_val.push_back(std::stod(s));
    }
    setting = new_val;
}