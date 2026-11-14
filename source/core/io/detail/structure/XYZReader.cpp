// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/detail/structure/XYZReader.h>

#include <io/File.h>
#include <settings/MoleculeSettings.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

io::pdb::PDBStructure io::detail::xyz::read(const io::File& path) {
    console::print_info("Reading XYZ structure file from \"" + path.str() + "\"");
    console::indent();
    io::pdb::PDBStructure res;
    if (settings::molecule::implicit_hydrogens) {
        console::print_text("Note: implicit hydrogens setting is enabled, but XYZ files do not contain enough information to assign hydrogens. Skipping them for this structure.");
    }
    res.supports_implicit_hydrogens = false;

    std::ifstream input(path);
    if (!input.is_open()) {throw except::io_error("XYZReader::read: Could not open file \"" + path.str() + "\"");}

    // header: the atom count, then a free-form comment line (in extended XYZ, the Properties and lattice keys)
    std::string line;
    int line_number = 1;
    auto error = [&] (const std::string& what) {
        return except::io_error("XYZReader::read: " + what + " on line " + std::to_string(line_number) + " of \"" + path.str() + "\"");
    };
    if (!std::getline(input, line)) {throw error("missing atom count");}
    auto count_token = utility::remove_leading_and_trailing(line, " \t\r");
    if (!utility::isinteger(count_token) || std::stoi(std::string(count_token)) <= 0) {
        throw error("expected a positive atom count, got \"" + std::string(count_token) + "\"");
    }
    int count = std::stoi(std::string(count_token));
    std::getline(input, line);
    ++line_number;

    // atoms: an element symbol and three coordinates, then any per-atom columns (tags, forces, ...), which are ignored
    res.atoms.reserve(count);
    while (static_cast<int>(res.atoms.size()) < count) {
        ++line_number;
        if (!std::getline(input, line)) {throw error("expected " + std::to_string(count) + " atoms, found " + std::to_string(res.atoms.size()));}
        auto tokens = utility::split(line, " \t\r");
        if (tokens.size() < 4 || !utility::isnumeric(tokens[1]) || !utility::isnumeric(tokens[2]) || !utility::isnumeric(tokens[3])) {
            throw error("expected an element symbol and three coordinates, got \"" + line + "\"");
        }

        if (!constants::symbols::detail::string_to_atomt_map.contains(tokens[0])) {throw error("unknown element \"" + tokens[0] + "\"");}
        PDBAtom atom;
        atom.set_element(constants::symbols::parse_element_string(tokens[0]));
        atom.coordinates().x() = std::stod(tokens[1]);
        atom.coordinates().y() = std::stod(tokens[2]);
        atom.coordinates().z() = std::stod(tokens[3]);
        res.atoms.emplace_back(std::move(atom));
    }

    // only single-frame files: anything but blank lines after the atoms is either another frame or a wrong atom count
    while (std::getline(input, line)) {
        ++line_number;
        if (!utility::remove_leading_and_trailing(line, " \t\r").empty()) {
            throw error("unexpected content after the " + std::to_string(count) + " declared atoms (multi-frame XYZ files are not supported)");
        }
    }

    console::print_text("Successfully read " + std::to_string(res.atoms.size()) + " atomic records.");
    console::unindent();
    return res;
}