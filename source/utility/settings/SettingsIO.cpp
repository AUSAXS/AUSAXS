#include <utility/settings/SettingsIO.h>
#include <utility/settings/SettingsIORegistry.h>
#include <utility/Utility.h>
#include <utility/Exceptions.h>

#include <fstream>
#include <filesystem>

void settings::detail::parse_option(const std::string& name, const std::vector<std::string>& value) {
    if (io::detail::settings_storage.count(name) == 0) {
        throw std::runtime_error("Unknown option: \"" + name + "\".");
    }
    io::detail::settings_storage[name]->set(value);
}

bool settings::detail::is_comment_char(char c) {
    switch (c) {
        case '#':
        case ';':
        case '/':
        case '[':
            return true;
        default:
            return false;
    }
}

void settings::read(const std::string& path) {
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("Settings::read: Could not open setup file.");}

    std::string line; 
    while (getline(input, line)) {
        if (line.empty()) {continue;}             // skip empty lines
        if (detail::is_comment_char(line[0])) {continue;} // skip comments

        auto tokens = utility::split(line, " \t");
        if (tokens.size() == 1) {throw except::io_error("settings::read: Invalid line in setup file \"" + line + "\": \"" + tokens[0] + "\"");}
        std::string name = tokens[0];
        tokens.erase(tokens.begin());
        detail::parse_option(name, tokens);
    }
}

void settings::write(const std::string& path) {
    utility::create_directory(path);
    std::ofstream output(path);
    if (!output.is_open()) {throw std::ios_base::failure("Settings::read: Could not open setup file.");}

    output << "### Auto-generated settings file ###\n";
    for (const auto& section : io::detail::settings_sections) {
        output << "\n[" << section.name << "]\n";
        for (const auto& setting : section.settings) {
            output << setting->names.front() << " " << setting->get() << std::endl;
        }
    }
}

bool settings::discover(std::string path) {
    static std::vector<std::string> valid_names = {"settings", "setting", "setup", "config"};
    if (path.back() != '/') {path += "/";}
    for (const auto& e : valid_names) {
        if (std::filesystem::exists(path + e + ".txt")) {
            settings::read(path + e + ".txt");
            return true;
        }
    }
    return false;
}