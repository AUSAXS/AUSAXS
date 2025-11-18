// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/detail/structure/XYZReader.h>
#include <io/ExistingFile.h>
#include <utility/Console.h>
#include <settings/MoleculeSettings.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

io::pdb::PDBStructure io::detail::xyz::read(const io::File& path) {
    console::print_info("Reading XYZ structure file from \"" + path.str() + "\"");
    console::indent();
    if (settings::molecule::implicit_hydrogens) {
        console::print_text("Note: implicit hydrogens setting is enabled, but XYZ files do not contain enough information to assign hydrogens. Disabling.");
        settings::molecule::implicit_hydrogens = false;
    }

    io::pdb::PDBStructure res;
    std::ifstream input(path);
    if (!input.is_open()) {throw except::io_error("XYZReader::read: Could not open file \"" + path.str() + "\"");}

    int section_id = 0; // 0 = header, 1 = atoms, 2 = footer
    std::string line;
    while(getline(input, line)) {
        if (utility::remove_all(line, " \n\r").empty()) {continue;}
        auto tokens = utility::split(line, " \t");

        // skip header
        if (tokens.size() != 4 && section_id == 0) {continue;}
        else {section_id = 1;}                              // progress to main data section
        if (tokens.size() != 4 && section_id == 1) {break;} // progress to footer

        PDBAtom atom;
        atom.element = constants::symbols::parse_element_string(tokens[0]);
        atom.coordinates().x() = std::stod(tokens[1]);
        atom.coordinates().y() = std::stod(tokens[2]);
        atom.coordinates().z() = std::stod(tokens[3]);
        res.atoms.emplace_back(std::move(atom));
    }

    console::print_text("Successfully read " + std::to_string(res.atoms.size()) + " atomic records.");
    console::unindent();
    return res;
}