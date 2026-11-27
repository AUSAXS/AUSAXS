// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <residue/detail/Residue.h>

#include <constants/Constants.h>
#include <io/ExistingFile.h>
#include <io/detail/structure/CIFReader.h>
#include <residue/detail/ResidueMap.h>
#include <utility/StringUtils.h>

#include <sstream>
#include <utility>

using namespace ausaxs;
using namespace ausaxs::residue::detail;

//### ATOM ###//
ResidueAtom::ResidueAtom(std::string name, std::string altname, constants::atom_t symbol) 
    : name(std::move(name)), altname(std::move(altname)), atom(symbol)
{}

void ResidueAtom::add_bond(const constants::atom_t symbol) {
    if (constants::symbols::is_hydrogen(symbol)) {
        hydrogen_bonds++;
    }
}

std::string ResidueAtom::to_string() const {
    return "Atom " + name + " " + altname + " with " + std::to_string(hydrogen_bonds) + " hydrogen bonds";
}


//### BOND ###//
Bond::Bond(std::string  name1, std::string  name2, int order) : name1(std::move(name1)), name2(std::move(name2)), order(order) {}

std::string Bond::to_string() const {
    return "Bond " + name1 + (order == 1 ? " - " : " = ") + name2;
}

int Bond::parse_order(const std::string& order) {
    std::string order_lower = utility::to_lowercase(order);
    if (order_lower == "sing") {return 1;}
    if (order_lower == "doub") {return 2;}
    if (order_lower == "trip") {return 3;}
    throw ausaxs::except::runtime_error("Bond::parse_order: Invalid bond order: " + order);
}

//### RESIDUE ###//
Residue::Residue(std::string  name) : name(std::move(name)) {}

Residue::Residue(std::string  name, std::vector<ResidueAtom> atoms, const std::vector<Bond>& bonds) : name(std::move(name)), atoms(std::move(atoms)) {
    apply_bond(bonds);
}

void Residue::add_atom(const std::string& name, const std::string& altname, constants::atom_t atom) {
    name_map.insert({name, atoms.size()});
    atoms.emplace_back(name, altname, atom);
}

void Residue::apply_bond(const std::vector<Bond>& bonds) {
    for (const Bond& b : bonds) {
        apply_bond(b);
    }
}

void Residue::apply_bond(const Bond& bond) {
    ResidueAtom& a1 = atoms.at(name_map.at(bond.name1));
    ResidueAtom& a2 = atoms.at(name_map.at(bond.name2));
    a1.add_bond(a2.atom);
    a2.add_bond(a1.atom);
}

std::string Residue::to_string() const {
    std::stringstream ss;        
    for (const ResidueAtom& a : atoms) {
        ss << a << std::endl;
    }
    return ss.str();
}

ResidueMap Residue::to_map() const {
    ResidueMap map;
    for (const ResidueAtom& a : atoms) {
        // skip all H's, they are automatically handled by the SimpleResidueMap
        if (constants::symbols::is_hydrogen(a.atom)) {
            continue;
        }

        // check if the alternate name should also be inserted
        if (a.altname != a.name && !a.altname.empty()) {
            map.insert(a.altname, a.atom, a.hydrogen_bonds);
        }
        map.insert(a.name, a.atom, a.hydrogen_bonds);
    }
    return map;
}

Residue Residue::parse(const io::ExistingFile& filename) {
    auto residues = io::detail::cif::read_residue(filename);
    if (1 < residues.size()) {throw except::io_error("Residue::parse: Expected a single residue in file \"" + filename.str() + "\"");}
    return residues.front();
}

//### STREAM OPERATORS ###//
std::ostream& residue::detail::operator<<(std::ostream& os, const ResidueAtom& a) {os << a.to_string(); return os;}
std::ostream& residue::detail::operator<<(std::ostream& os, const Bond& b) {os << b.to_string(); return os;}
std::ostream& residue::detail::operator<<(std::ostream& os, const Residue& l) {os << l.to_string(); return os;}