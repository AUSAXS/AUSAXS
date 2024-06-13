/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/gmx/SmartOption.h>

#include <fstream>
#include <sstream>
#include <filesystem>

template<> void md::SmartOption<double>::set(const std::string& value) {this->value = std::stod(value);}
template<> void md::SmartOption<int>::set(const std::string& value) {this->value = std::stoi(value);}
template<> void md::SmartOption<std::string>::set(const std::string& value) {this->value = value;}
template<> std::string md::SmartOption<double>::get() const {return std::to_string(value);}
template<> std::string md::SmartOption<int>::get() const {return std::to_string(value);}
template<> std::string md::SmartOption<std::string>::get() const {return value;}

std::unordered_map<std::string, md::ISmartOption*>& md::ISmartOption::get_all_options() {
    static std::unordered_map<std::string, ISmartOption*> all_options;
    return all_options;
}

void md::settings::parse_option(const std::string& name, const std::string& value) {
    if (ISmartOption::get_all_options().count(name) == 0) {
        throw std::runtime_error("Unknown option: " + name);
    }
    ISmartOption::get_all_options()[name]->set(value);
}

bool md::settings::is_comment_char(char c) {
    switch (c) {
        case '#':
        case ';':
        case '/':
            return true;
        default:
            return false;
    }
}

void md::settings::read(const std::string& filename) {
    std::ifstream in(filename);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || is_comment_char(line[0])) {continue;}
        std::string name, value;
        std::istringstream iss(line);
        iss >> name >> value;
        parse_option(name, value);
    }
}

void md::settings::write(const std::string& filename) {
    std::ofstream out(filename);
    for (const auto& [name, option] : ISmartOption::get_all_options()) {
        out << name << " " << option->get() << std::endl;
    }
}

bool md::settings::discover(std::string path) {
    if (path[path.size()-1] != '/') {path += "/";}
    std::vector<std::string> valid_names = {"settings", "setting", "setup", "config"};
    for (const auto& e : valid_names) {
        if (std::filesystem::exists(path + e + ".txt")) {
            std::cout << "Using discovered settings file at " << path + e + ".txt" << std::endl;
            settings::read(path + e + ".txt");
            return true;
        }
    }
    return false;
}
