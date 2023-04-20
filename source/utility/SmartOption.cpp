#include <utility/SmartOption.h>
#include <utility/Exceptions.h>

#include <fstream>
#include <sstream>
#include <filesystem>
#include <iostream>

template<> std::string settings::detail::SmartOption<std::string>::get() const {return value;}
template<> std::string settings::detail::SmartOption<double>::get() const {return std::to_string(value);}
template<> std::string settings::detail::SmartOption<int>::get() const {return std::to_string(value);}
template<> std::string settings::detail::SmartOption<unsigned int>::get() const {return std::to_string(value);}
template<> std::string settings::detail::SmartOption<bool>::get() const {return std::to_string(value);}

template<> std::string settings::detail::SmartOption<std::vector<std::string>>::get() const {
    std::string str;
    std::for_each(value.begin(), value.end(), [&str] (const std::string& s) {str += s + " ";});
    return str;
}

template<> std::string settings::detail::SmartOption<std::vector<double>>::get() const {
    std::string str;
    std::for_each(value.begin(), value.end(), [&str] (double s) {str += std::to_string(s) + " ";});
    return str;
}

template<> std::string settings::detail::SmartOption<std::vector<int>>::get() const {
    std::string str;
    std::for_each(value.begin(), value.end(), [&str] (int s) {str += std::to_string(s) + " ";});
    return str;
}

template<> void settings::detail::SmartOption<std::string>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many values.");}
    value = str[0];
}

template<> void settings::detail::SmartOption<bool>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many values.");}

    if (str[0] == "true" || str[0] == "TRUE" || str[0] == "1") {value = true; return;}
    else if (str[0] == "false" || str[0] == "FALSE" || str[0] == "0") {value = false; return;}
    throw except::parse_error("Settings::parse_bool: Option \"" + get() + "\" expected boolean string, but got \"" + str[0] + "\".");
}

template<> void settings::detail::SmartOption<double>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many values.");}
    value = std::stod(str[0]);
}

template<> void settings::detail::SmartOption<int>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many values.");}
    value = std::stoi(str[0]); 
}

template<> void settings::detail::SmartOption<unsigned int>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many values.");}
    value = std::stoi(str[0]);
}

template<> void settings::detail::SmartOption<std::vector<std::string>>::set(const std::vector<std::string>& str) {
    value = str;
}

template<> void settings::detail::SmartOption<std::vector<double>>::set(const std::vector<std::string>& str) {
    std::vector<double> new_val;
    for (auto& s : str) {
        if (s.empty() || s == " ") {continue;}
        new_val.push_back(std::stod(s));
    }
    if (new_val.empty()) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received no values.");}
    value = new_val;
}

template<> void settings::detail::SmartOption<std::vector<int>>::set(const std::vector<std::string>& str) {
    std::vector<int> new_val;
    for (auto& s : str) {
        if (s.empty() || s == " ") {continue;}
        new_val.push_back(std::stoi(s));
    }
    if (new_val.empty()) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received no values.");}
    value = new_val;
}