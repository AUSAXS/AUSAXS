// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <residue/ResidueStorage.h>
#include <residue/detail/ResidueMap.h>
#include <residue/detail/InvalidResidueMap.h>
#include <residue/detail/ResidueStorageBasis.h>
#include <io/ExistingFile.h>
#include <utility/Curl.h>
#include <utility/Console.h>
#include <settings/GeneralSettings.h>
#include <constants/Constants.h>

#include <fstream>
#include <filesystem>
#include <unordered_map>
#include <regex>
#include <iostream>

using namespace ausaxs;
using namespace ausaxs::residue;
using namespace ausaxs::residue::detail;

ResidueStorage::ResidueStorage() = default;

ResidueStorage::~ResidueStorage() = default;

void ResidueStorage::insert(const std::string& name, const ResidueMap& residue) {
    data.emplace(name, residue);
}

bool ResidueStorage::contains(const std::string& name) {
    if (!initialized) {initialize();}
    return data.contains(name);
}

ResidueMap& ResidueStorage::get(const std::string& name) {
    if (!initialized) {initialize();}
    if (data.find(name) == data.end()) {
        bool downloaded = false;

        // small cache to avoid spamming the console with the same download
        static std::unordered_map<std::string, bool> seen_before;
        if (!seen_before.contains(name)) {
            console::indent(2);
            console::print_text("Unknown residue: \"" + name + "\". Attempting to download specification.");
            console::indent();
            seen_before[name] = true;
            downloaded = update_or_download_residue(name);
            console::unindent(3);
        }

        if (!downloaded) {
            // singleton for all unknown residues which will return 0 hydrogens for all atoms
            static InvalidResidueMap invalid_residue;
            return invalid_residue;
        }
    }
    return data.at(name);
}

void ResidueStorage::initialize() {
    initialized = true;
    std::string path = settings::general::residue_folder;
    io::Folder(path).create();
    std::ifstream file(path + "master.dat");
    if (!file.is_open()) {
        write_master_basis();
        file.open(path + "master.dat");
        if (!file.is_open()) {throw except::io_error("ResidueStorage::initialize: Could not open file: " + path + "master" + ".dat");}
    }

    std::string line;
    while (file.peek() != EOF) {
        std::getline(file, line);

        // skip until we reach the start of a residue
        if (line.find("#") == std::string::npos) {
            continue;
        } else {
            // the line following the # is the name of the residue
            std::getline(file, line);
            std::string residue = line;

            // prepare map
            std::unordered_map<AtomKey, int> map;
            while (file.peek() != EOF) {
                std::getline(file, line);
                // stop if we reach the start of a new residue
                if (line.empty() || line.find("#") != std::string::npos) {
                    break;
                }

                // lines are of the form "element atom hydrogens"
                std::vector<std::string> tokens = utility::split(line, " \n\r");
                if (tokens.size() != 3) {throw except::io_error("ResidueStorage::initialize: Invalid line in master file: " + line + ". Perhaps the file is corrupted. Delete it to regenerate.");}
                std::string element = tokens[0];
                std::string atom = tokens[1];
                int hydrogens = std::stoi(tokens[2]);
                map.emplace(AtomKey(atom, constants::symbols::parse_element_string(element)), hydrogens);
            }
            insert(residue, std::move(map));
        }
    }
}

bool ResidueStorage::update_or_download_residue(const std::string& name) {
    std::string path = settings::general::residue_folder;
    std::regex regex("[A-Z0-9]{1,3}");

    if (std::regex_match(name, regex)) {
        // check if the file already exists. if not, download it.
        if (!std::filesystem::exists(path + name + ".cif")) {
            if (settings::general::offline) {
                console::print_warning("ResidueStorage::download_residue: \"" + name + "\" specification cannot be downloaded as offline mode is enabled!");
                return false;
            }

            bool downloaded = curl::download("files.rcsb.org/ligands/view/" + name + ".cif", path + name + ".cif"); // download the cif file
            if (!downloaded) {return false;}
        } else {
            console::print_text_minor("Residue " + name + " is already downloaded, but not present in the master list.");
            console::print_text_minor("The file will now be parsed and re-added to the master file.");
        }

        // parse the cif file & add it to storage
        residue::detail::ResidueMap map;
        try {
            map = Residue::parse(path + name + ".cif").to_map();
        } catch (const std::exception& e) {
            // if the residue could not be parsed, try to redownload it and parse it again
            console::print_warning("ResidueStorage::download_residue: Could not parse residue: " + name + ". \nError: " + e.what());
            console::print_text("Attempting to redownload the residue definition.");
            curl::download("files.rcsb.org/ligands/view/" + name + ".cif", path + name + ".cif");
            map = Residue::parse(path + name + ".cif").to_map();
        }
        insert(name, map);

        // write the residue to the master file
        write_residue(name);
    } else {
        throw except::io_error("ResidueStorage::download_residue: Invalid residue name: \"" + name + "\"; expected a 2-3 letter code.");
    }
    return true;
}

void ResidueStorage::write_residue(const std::string& name) {
    std::string path = settings::general::residue_folder;
    io::Folder(path).create();
    std::ofstream file(path + "master.dat", std::ios::app); // open in append mode
    if (!file.is_open()) {throw except::io_error("ResidueStorage::write_residue: Could not open file: " + path + "master" + ".dat");}

    // write the map to the master file
    auto map = get(name);
    file << "#" << "\n" << name << "\n"; // residue header
    for (const auto& [key, val] : map.get_backing_map()) {
        file << constants::symbols::to_string(key.atom) << " " << key.name << " " << val << "\n";
    }
    file << std::endl;
}

constants::atomic_group_t ResidueStorage::get_atomic_group(const std::string& residue_name, const std::string& atom_name, constants::atom_t atom) {
    if (!initialized) {initialize();}
    auto& residue = get(residue_name);
    return residue.get_atomic_group(atom_name, atom);
}