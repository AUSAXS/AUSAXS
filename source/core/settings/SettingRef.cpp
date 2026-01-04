// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/SettingRef.h>
#include <utility/Exceptions.h>
#include <utility/Limit.h>

#include <algorithm>

using namespace ausaxs;

std::unordered_map<std::string, std::shared_ptr<settings::io::detail::ISettingRef>>& settings::io::detail::ISettingRef::get_stored_settings() {
    static std::unordered_map<std::string, std::shared_ptr<settings::io::detail::ISettingRef>> stored_settings;
    return stored_settings;
}

template<> std::string ausaxs::settings::io::detail::type_as_string<std::string>(const std::string&) {return "string";}
template<> std::string ausaxs::settings::io::detail::type_as_string<double>(const double&) {return "double";}
template<> std::string ausaxs::settings::io::detail::type_as_string<int>(const int&) {return "int";}
template<> std::string ausaxs::settings::io::detail::type_as_string<unsigned int>(const unsigned int&) {return "uint";}
template<> std::string ausaxs::settings::io::detail::type_as_string<bool>(const bool&) {return "bool";}
template<> std::string ausaxs::settings::io::detail::type_as_string<std::vector<std::string>>(const std::vector<std::string>&) {return "vector-string";}
template<> std::string ausaxs::settings::io::detail::type_as_string<std::vector<double>>(const std::vector<double>&) {return "vector-double";}
template<> std::string ausaxs::settings::io::detail::type_as_string<std::vector<int>>(const std::vector<int>&) {return "vector-int";}
template<> std::string ausaxs::settings::io::detail::type_as_string<ausaxs::Limit>(const ausaxs::Limit&) {return "limit";}

template<> std::string settings::io::detail::SettingRef<std::string>::get() const {return settingref;}
template<> std::string settings::io::detail::SettingRef<double>::get() const {return std::to_string(settingref);}
template<> std::string settings::io::detail::SettingRef<int>::get() const {return std::to_string(settingref);}
template<> std::string settings::io::detail::SettingRef<unsigned int>::get() const {return std::to_string(settingref);}
template<> std::string settings::io::detail::SettingRef<bool>::get() const {return std::to_string(settingref);}
template<> std::string settings::io::detail::SettingRef<Limit>::get() const {return std::to_string(settingref.min) + " " + std::to_string(settingref.max);}
template<> std::string settings::io::detail::SettingRef<std::vector<std::string>>::get() const {
    std::string str;
    std::for_each(settingref.begin(), settingref.end(), [&str] (const std::string& s) {str += s + " ";});
    return str;
}
template<> std::string settings::io::detail::SettingRef<std::vector<double>>::get() const {
    std::string str;
    std::for_each(settingref.begin(), settingref.end(), [&str] (double s) {str += std::to_string(s) + " ";});
    return str;
}
template<> std::string settings::io::detail::SettingRef<std::vector<int>>::get() const {
    std::string str;
    std::for_each(settingref.begin(), settingref.end(), [&str] (int s) {str += std::to_string(s) + " ";});
    return str;
}


template<> void settings::io::detail::SettingRef<std::string>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many settings.");}
    settingref = str[0];
}
template<> void settings::io::detail::SettingRef<bool>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many settings.");}

    if (str[0] == "true" || str[0] == "TRUE" || str[0] == "1") {settingref = true; return;}
    else if (str[0] == "false" || str[0] == "FALSE" || str[0] == "0") {settingref = false; return;}
    throw except::parse_error("Settings::parse_bool: Option \"" + get() + "\" expected boolean string, but got \"" + str[0] + "\".");
}
template<> void settings::io::detail::SettingRef<double>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many settings.");}
    settingref = std::stod(str[0]);
}
template<> void settings::io::detail::SettingRef<int>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many settings.");}
    settingref = std::stoi(str[0]); 
}
template<> void settings::io::detail::SettingRef<unsigned int>::set(const std::vector<std::string>& str) {
    if (str.size() != 1) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many settings.");}
    settingref = std::stoi(str[0]);
}
template<> void settings::io::detail::SettingRef<std::vector<std::string>>::set(const std::vector<std::string>& str) {
    settingref = str;
}

template<> void settings::io::detail::SettingRef<Limit>::set(const std::vector<std::string>& str) {
    if (str.size() != 2) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received too many settings.");}
    settingref.min = std::stod(str[0]);
    settingref.max = std::stod(str[1]);
}

template<> void settings::io::detail::SettingRef<std::vector<double>>::set(const std::vector<std::string>& str) {
    std::vector<double> new_val;
    for (auto& s : str) {
        if (s.empty() || s == " ") {continue;}
        new_val.push_back(std::stod(s));
    }
    if (new_val.empty()) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received no settings.");}
    settingref = std::move(new_val);
}
template<> void settings::io::detail::SettingRef<std::vector<int>>::set(const std::vector<std::string>& str) {
    std::vector<int> new_val;
    for (auto& s : str) {
        if (s.empty() || s == " ") {continue;}
        new_val.push_back(std::stoi(s));
    }
    if (new_val.empty()) {throw except::parse_error("Settings::SmartOption::parse: Option \"" + get() + "\" received no settings.");}
    settingref = std::move(new_val);
}