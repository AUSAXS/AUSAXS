#include <residue/ResidueParser.h>
#include <residue/ResidueMap.h>
#include <utility/Exceptions.h>
#include <constants/Constants.h>
#include <io/ExistingFile.h>

#include <sstream>
#include <fstream>

using namespace residue::detail;

//### ATOM ###//
Atom::Atom(const std::string& name, const std::string& altname, constants::atom_t atom) : name(name), altname(altname), atom(atom) {
    valency = constants::valence::get_valence(atom);
}

Atom::Atom(const std::string& name, int charge, constants::atom_t atom) : name(name), atom(atom) {
    // the goal of this whole class is to determine the total charge surrounding an atom
    // we do this by counting the "hidden" hydrogen bonds not typically present in a PDB file
    // thus the number of hydrogen bonds is later used as the effective charge of the atom
    // and so, when dealing with an ion, we just set the hydrogen_bonds equal to the charge
    // this is a hack solution, but it works as it should
    hydrogen_bonds = charge; // adding additional hydrogen bonds is directly translated to adding additional charge
}

void Atom::add_bond(const constants::atom_t atom, unsigned int order) {
    if (atom == constants::atom_t::H) {
        hydrogen_bonds++;
    }
    valency -= order;
}

std::string Atom::to_string() const {
    return "Atom " + name + " with valency " + std::to_string(valency) + " and " + std::to_string(hydrogen_bonds) + " hydrogen bonds";
}


//### BOND ###//
Bond::Bond(const std::string& name1, const std::string& name2, unsigned int order) : name1(name1), name2(name2), order(order) {}

std::string Bond::to_string() const {
    return "Bond " + name1 + (order == 1 ? " - " : " = ") + name2;
}

unsigned int Bond::parse_order(const std::string& order) {
    if (order == "SING") {return 1;}
    else if (order == "DOUB") {return 2;}
    else if (order == "TRIP") {return 3;}
    else {throw std::runtime_error("Bond::parse_order: Invalid bond order: " + order);}
}


//### RESIDUE ###//
Residue::Residue(const std::string& name) : name(name) {}

Residue::Residue(const std::string& name, std::vector<Atom> atoms, std::vector<Bond> bonds) : name(name), atoms(atoms) {
    apply_bond(bonds);
}

void Residue::add_atom(const std::string& name, const std::string& altname, constants::atom_t atom) {
    name_map.insert({name, atoms.size()});
    atoms.push_back(Atom(name, altname, atom));
}

void Residue::add_atom(const std::string& name, int charge, constants::atom_t atom) {
    name_map.insert({name, atoms.size()});
    atoms.push_back(Atom(name, charge, atom));
}

void Residue::apply_bond(const std::vector<Bond>& bonds) {
    for (const Bond& b : bonds) {
        apply_bond(b);
    }
}

void Residue::apply_bond(const Bond& bond) {
    Atom& a1 = atoms.at(name_map.at(bond.name1));
    Atom& a2 = atoms.at(name_map.at(bond.name2));
    a1.add_bond(a2.atom, bond.order);
    a2.add_bond(a1.atom, bond.order);
}

std::string Residue::to_string() const {
    std::stringstream ss;        
    for (const Atom& a : atoms) {
        ss << a << std::endl;
    }
    return ss.str();
}

ResidueMap Residue::to_map() const {
    ResidueMap map;
    for (const Atom& a : atoms) {
        // skip all H's, they are automatically handled by the SimpleResidueMap
        if (a.atom == constants::atom_t::H) {
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
    std::ifstream file(filename);

    std::string line;
    Residue residue(filename.stem());
    bool found_atom_section = false, found_bond_section = false;

    // check if we are dealing with a single ion
    while (std::getline(file, line)) {
        if (line.find("formula") != std::string::npos) {
            std::vector<std::string> tokens = utility::split(line, " \n\r");
            std::string formula = tokens[1];

            // the formula is of the form "Xn Ym" for residues (e.g. "C10 H22 O6"), but X for ions (e.g. "Zn"). 
            // if it is an ion, X can be parsed as an element
            if (!constants::symbols::detail::string_to_atomt_map.contains(formula)) {
                break;
            }

            // parse ion
            while(std::getline(file, line)) {
                // find the formal charge
                if (line.find("pdbx_formal_charge") != std::string::npos) {
                    tokens = utility::split(line, " \n\r");
                    int charge = std::stoi(tokens[1]);
                    residue.add_atom(formula, charge, constants::symbols::parse_element_string(formula));
                    return residue;
                }
            }
        }
    }

    // find the beginning of the atom section
    while (std::getline(file, line)) {
        if (line.find("atom.pdbx_ordinal") != std::string::npos) {
            found_atom_section = true;
            break;
        }
    }

    // parse the atoms in the atom section
    while (std::getline(file, line)) {
        // check for end of section
        if (line.find("#") != std::string::npos) {
            break;
        }

        std::vector<std::string> tokens = utility::split(line, " \n\r");
        std::string atom_id = utility::remove_quotation_marks(tokens[1]);
        std::string atom_id_alt = utility::remove_quotation_marks(tokens[2]);
        std::string type_symbol = tokens[3];

        // Sometimes the "1" in e.g. CD1 is omitted
        if (atom_id == atom_id_alt) {
            if (auto N = atom_id.size(); N > 2) { // second character must be a locator e.g. A, B, C, D, E, ...
                if (atom_id[N - 1] == '1' && !std::isdigit(atom_id[N - 2])) { // last character must be 1 & previous character must be a locator
                    atom_id_alt = atom_id.substr(0, N - 1);
                }
            }
        }

        residue.add_atom(atom_id, atom_id_alt, constants::symbols::parse_element_string(type_symbol));
    }

    // find the beginning of the bond section
    while (std::getline(file, line)) {
        if (line.find("bond.pdbx_ordinal") != std::string::npos) {
            found_bond_section = true;
            break;
        }
    }

    // parse the bonds in the bond section
    while (std::getline(file, line)) {
        // check for end of section
        if (line.find("#") != std::string::npos) {
            break;
        }

        std::vector<std::string> tokens = utility::split(line, " \n\r");
        std::string atom1 = utility::remove_quotation_marks(tokens[1]);
        std::string atom2 = utility::remove_quotation_marks(tokens[2]);
        unsigned int order = Bond::parse_order(tokens[3]);
        residue.apply_bond(Bond(atom1, atom2, order));
    }

    // remove a hydrogen bond from the N-terminus since it will almost always be bonded to the CA of the next chain
    if (residue.name_map.contains("N")) {
        residue.atoms[(residue.name_map.at("N"))].hydrogen_bonds -= 1;
    }

    // check that the file was read correctly
    if (!found_atom_section || !found_bond_section) {
        throw except::io_error("Could not find atom or bond section in file \"" + filename + "\"");
    }
    
    return residue;
}

//### STREAM OPERATORS ###//
std::ostream& residue::detail::operator<<(std::ostream& os, const Atom& a) {os << a.to_string(); return os;}
std::ostream& residue::detail::operator<<(std::ostream& os, const Bond& b) {os << b.to_string(); return os;}
std::ostream& residue::detail::operator<<(std::ostream& os, const Residue& l) {os << l.to_string(); return os;}