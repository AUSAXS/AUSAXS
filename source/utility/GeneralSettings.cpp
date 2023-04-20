#include <utility/GeneralSettings.h>
#include <utility/Utility.h>
#include <utility/Exceptions.h>

#include <fstream>
#include <filesystem>
#include <iostream>
#include <thread>

using namespace settings::detail;

namespace settings::general {
    SmartOption<std::string> residue_folder("temp/residues/", "residue-folder");
    SmartOption<bool> verbose(true, "verbose");
    SmartOption<unsigned int> threads(std::thread::hardware_concurrency(), "threads");
    SmartOption<std::string> output("output/", {"output", "output-folder"});
    SmartOption<bool> keep_hydrogens(false, "keep-hydrogens");
    SmartOption<bool> supplementary_plots(true, {"extra-plots", "additional-plots"});

    namespace detail {
        SmartOption<unsigned int> job_size(200); // The number of atoms to process in each job.
    };
}

void settings::detail::parse_option(const std::string& name, const std::vector<std::string>& value) {
    if (detail::ISmartOption::all_options.count(name) == 0) {
        throw std::runtime_error("Unknown option: " + name);
    }
    detail::ISmartOption::all_options[name]->set(value);
}

bool settings::detail::is_comment_char(char c) {
    switch (c) {
        case '#':
        case ';':
        case '/':
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
        if (is_comment_char(line[0])) {continue;} // skip comments

        auto tokens = utility::split(line, " \t");
        if (tokens.size() == 1) {throw except::io_error("settings::read: Invalid line in setup file \"" + line + "\": \"" + tokens[0] + "\"");}
        std::string name = tokens[0];
        tokens.erase(tokens.begin());
        parse_option(name, tokens);
    }
}

void settings::write(const std::string& path) {
    utility::create_directory(path);
    std::ofstream output(path);
    if (!output.is_open()) {throw std::ios_base::failure("Settings::read: Could not open setup file.");}

    output << "### Auto-generated settings file ###\n";
    for (const auto& [name, option] : detail::ISmartOption::all_options) {
        output << name << " " << option->get() << std::endl;
    }
}

bool settings::discover(std::string path) {
    static std::vector<std::string> valid_names = {"settings", "setting", "setup", "config"};
    if (path.back() != '/') {path += "/";}
    for (const auto& e : valid_names) {
        if (std::filesystem::exists(path + e + ".txt")) {
            std::cout << "Using discovered settings file at " << path + e + ".txt" << std::endl;
            settings::read(path + e + ".txt");
            return true;
        }
    }
    return false;
}