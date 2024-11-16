/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <residue/ResidueStorage.h>
#include <residue/ResidueMap.h>
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

ResidueStorage::ResidueStorage() {}

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
        console::print_info("Unknown residue: \"" + name + "\". Attempting to download specification.");
        download_residue(name);
    }
    return data.at(name);
}

void ResidueStorage::initialize() {
    initialized = true;
    std::string path = settings::general::residue_folder;
    io::Folder(path).create();
    std::ifstream file(path + "master.dat");
    if (!file.is_open()) {
        return;
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

void ResidueStorage::download_residue(const std::string& name) {
    std::string path = settings::general::residue_folder;
    std::regex regex("[A-Z0-9]{2,3}");

    if (std::regex_match(name, regex)) {
        // check if the file already exists. if not, download it.
        if (!std::filesystem::exists(path + name + ".cif")) {
            curl::download("files.rcsb.org/ligands/view/" + name + ".cif", path + name + ".cif"); // download the cif file
        } else {
            std::cout << "\tResidue " << name << " is already downloaded, but not present in the master list. \n\tThe file will now be parsed and re-added to the master file." << std::endl;
        }

        // parse the cif file & add it to storage
        insert(name, Residue::parse(path + name + ".cif").to_map());

        // write the residue to the master file
        write_residue(name);
    } else {
        throw except::io_error("ResidueStorage::download_residue: Invalid residue name: \"" + name + "\"; expected a 2-3 letter code.");
    }
}

void ResidueStorage::write_residue(const std::string& name) {
    std::string path = settings::general::residue_folder;
    io::Folder(path).create();
    std::ofstream file(path + "master.dat", std::ios::app); // open in append mode
    if (!file.is_open()) {throw except::io_error("ResidueStorage::write_residue: Could not open file: " + path + "master" + ".dat");}

    // write the map to the master file
    auto map = get(name);
    file << "#" << "\n" << name << "\n"; // residue header
    for (const auto& [key, val] : map) {
        file << constants::symbols::to_string(key.atom) << " " << key.name << " " << val << "\n";
    }
    file << std::endl;
}

constants::atomic_group_t ResidueStorage::get_atomic_group(const std::string& residue_name, const std::string& atom_name, constants::atom_t atom) {
    if (!initialized) {initialize();}
    auto& residue = get(residue_name);
    return residue.get_atomic_group(atom_name, atom);
}