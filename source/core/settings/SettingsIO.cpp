// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/SettingsIO.h>
#include <settings/SettingsIORegistry.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <utility/Console.h>

#include <fstream>

using namespace ausaxs;

void settings::detail::parse_option(const std::string& name, const std::vector<std::string>& value) {
    if (!settings::io::detail::ISettingRef::get_stored_settings().contains(name)) {
        throw std::runtime_error("Unknown option: \"" + name + "\".");
    }
    settings::io::detail::ISettingRef::get_stored_settings()[name]->set(value);
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

void settings::read(const ::io::ExistingFile& path) {
    console::print_info("Reading settings from file: \"" + path.str() + "\"");

    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("settings::read: Could not open setup file.");}

    std::string line; 
    while (getline(input, line)) {
        if (line.empty()) {continue;}                       // skip empty lines
        if (detail::is_comment_char(line[0])) {continue;}   // skip comments

        auto tokens = utility::split(line, " \t");
        if (tokens.size() == 1) {
            console::print_text_minor("\tsettings::read: Skipping line: \"" + line + "\" as no value could be found.");
            continue;
        }
        std::string name = tokens[0];
        tokens.erase(tokens.begin());
        detail::parse_option(name, tokens);
    }
}

void settings::write(const ::io::File& path) {
    path.directory().create();
    std::ofstream output(path);
    if (!output.is_open()) {throw std::ios_base::failure("settings::write: Could not open setup file.");}

    output << "### Auto-generated settings file ###\n";
    for (const auto& section : settings::io::SettingSection::get_sections()) {
        output << "\n[   " << section->name << "   ]\n";
        for (const auto& setting : section->settings) {
            output << setting->names.front() << " " << setting->get() << std::endl;
        }
    }
}

bool settings::discover(const ::io::Folder& path) {
    static std::vector<std::string> valid_names = {"settings", "setting", "setup", "config"};
    for (const auto& e : valid_names) {
        ::io::File file(path, e, ".txt");
        if (file.exists()) {
            settings::read(file);
            return true;
        }
    }
    return false;
}