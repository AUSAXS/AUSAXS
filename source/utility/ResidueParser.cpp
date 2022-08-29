#include <utility/ResidueParser.h>
#include <utility/Curl.h>
#include <utility/Settings.h>
#include <utility/Constants.h>

#include <regex>
#include <filesystem>
#include <fstream>

//### ATOM ###//
parser::residue::detail::Atom::Atom(std::string name, std::string altname, std::string symbol) : name(name), altname(altname), symbol(symbol) {
    valency = constants::valence::atomic.get(symbol);
}
                
parser::residue::detail::Atom::Atom(std::string name, int charge) : name(name) {
    // the goal of this whole class is to determine the total charge surrounding an atom
    // we do this by counting the "hidden" hydrogen bonds not typically present in a PDB file
    // thus the number of hydrogen bonds is later used as the effective charge of the atom
    // and so, when dealing with an ion, we just set the hydrogen_bonds equal to the charge
    // this is a hack solution, but it works as it should
    hydrogen_bonds = charge;
}

void parser::residue::detail::Atom::add_bond(std::string symbol, unsigned int order) {
    if (symbol == "H") {
        hydrogen_bonds++;
    }

    valency -= order;
}

std::string parser::residue::detail::Atom::to_string() const {
    return "Atom " + name + " with valency " + std::to_string(valency) + " and " + std::to_string(hydrogen_bonds) + " hydrogen bonds";
}


//### BOND ###//
parser::residue::detail::Bond::Bond(std::string name1, std::string name2, unsigned int order) : name1(name1), name2(name2), order(order) {}

std::string parser::residue::detail::Bond::to_string() const {
    return "Bond " + name1 + (order == 1 ? " - " : " = ") + name2;
}

unsigned int parser::residue::detail::Bond::parse_order(std::string order) {
    if (order == "SING") {return 1;}
    else if (order == "DOUB") {return 2;}
    else if (order == "TRIP") {return 3;}
    else {throw std::runtime_error("Invalid bond order: " + order);}
}


//### RESIDUE ###//
parser::residue::detail::Residue::Residue(std::string name) : name(name) {}

parser::residue::detail::Residue::Residue(std::string name, std::vector<Atom> atoms, std::vector<Bond> bonds) : name(name), atoms(atoms) {
    apply_bond(bonds);
}

void parser::residue::detail::Residue::add_atom(std::string name, std::string altname, std::string symbol) {
    name_map.insert({name, atoms.size()});
    atoms.push_back(Atom(name, altname, symbol));
}

void parser::residue::detail::Residue::add_atom(std::string name, int charge) {
    name_map.insert({name, atoms.size()});
    atoms.push_back(Atom(name, charge));
}

void parser::residue::detail::Residue::apply_bond(const std::vector<Bond>& bonds) {
    for (const Bond& b : bonds) {
        apply_bond(b);
    }
}

void parser::residue::detail::Residue::apply_bond(const Bond& bond) {
    Atom& a1 = atoms.at(name_map.at(bond.name1));
    Atom& a2 = atoms.at(name_map.at(bond.name2));
    a1.add_bond(a2.symbol, bond.order);
    a2.add_bond(a1.symbol, bond.order);
}

std::string parser::residue::detail::Residue::to_string() const {
    std::stringstream ss;        
    for (const Atom& a : atoms) {
        ss << a << std::endl;
    }
    return ss.str();
}

saxs::detail::SimpleResidueMap parser::residue::detail::Residue::to_map() const {
    saxs::detail::SimpleResidueMap map;
    for (const Atom& a : atoms) {
        // skip all H's, they are automatically handled by the SimpleResidueMap
        if (a.name.find("H") != std::string::npos) {
            continue;
        }

        // check if the alternate name should also be inserted
        if (a.altname != a.name && !a.altname.empty()) {
            map.insert(a.altname, a.hydrogen_bonds);
        }
        map.insert(a.name, a.hydrogen_bonds);
    }
    return map;
}

parser::residue::detail::Residue parser::residue::detail::Residue::parse(std::string filename) {
    std::ifstream file(filename);

    std::string line;
    Residue residue(utility::stem(filename));
    bool found_atom_section = false, found_bond_section = false;

    // check if we are dealing with a single ion
    while (std::getline(file, line)) {
        if (line.find("formula") != std::string::npos) {
            std::vector<std::string> tokens = utility::split(line, " \n\r");
            std::string formula = tokens[1];

            // the formula is of the form "Xn Ym" for residues, but X for ions. 
            // if it is an ion, X should be present in the mass map so we use this as a check
            if (!constants::mass::atomic.contains(formula)) {
                break;
            }

            // parse ion
            while(std::getline(file, line)) {
                // find the formal charge
                if (line.find("pdbx_formal_charge") != std::string::npos) {
                    tokens = utility::split(line, " \n\r");
                    int charge = std::stoi(tokens[1]);
                    residue.add_atom(formula, charge);

                    std::cout << "Parsed a single ion " << formula << " with charge " << charge << std::endl;
                    return residue;
                }
            }
        }
    }

    // find the beginning of the atom section
    while (std::getline(file, line)) {
        if (line.find("atom.pdbx_ordinal") != std::string::npos) {
            // std::cout << "Found start of atom section" << std::endl;
            found_atom_section = true;
            break;
        }
    }

    // parse the atoms in the atom section
    while (std::getline(file, line)) {
        // check for end of section
        if (line.find("#") != std::string::npos) {
            // std::cout << "Found end of atom section" << std::endl;
            break;
        }

        std::vector<std::string> tokens = utility::split(line, " \n\r");
        std::string atom_id = utility::remove_quotation_marks(tokens[1]);
        std::string atom_id_alt = utility::remove_quotation_marks(tokens[2]);
        std::string type_symbol = tokens[3];

        // HN is sometimes used as an alias for "H"
        // if (atom_id_alt == "H") {atom_id_alt = "HN";}

        // Sometimes the "1" in e.g. CD1 is omitted
        if (atom_id == atom_id_alt) {
            if (atom_id[atom_id.size() - 1] == '1') {
                atom_id_alt = atom_id.substr(0, atom_id.size()-1);
            }
        }

        residue.add_atom(atom_id, atom_id_alt, type_symbol);
    }

    // find the beginning of the bond section
    while (std::getline(file, line)) {
        if (line.find("bond.pdbx_ordinal") != std::string::npos) {
            // std::cout << "Found start of bond section" << std::endl;
            found_bond_section = true;
            break;
        }
    }

    // parse the bonds in the bond section
    while (std::getline(file, line)) {
        // check for end of section
        if (line.find("#") != std::string::npos) {
            // std::cout << "Found end of bond section" << std::endl;
            break;
        }

        std::vector<std::string> tokens = utility::split(line, " \n\r");
        std::string atom1 = utility::remove_quotation_marks(tokens[1]);
        std::string atom2 = utility::remove_quotation_marks(tokens[2]);
        unsigned int order = Bond::parse_order(tokens[3]);
        residue.apply_bond(Bond(atom1, atom2, order));
    }

    // check that the file was read correctly
    if (!found_atom_section || !found_bond_section) {
        throw except::io_error("Could not find atom or bond section in file " + filename);
    }
    
    return residue;
}

//### RESIDUE STORAGE ###//
parser::residue::ResidueStorage::ResidueStorage() {
    initialize();
}

void parser::residue::ResidueStorage::insert(std::string name, saxs::detail::SimpleResidueMap residue) {
    data.emplace(name, residue);
}

saxs::detail::SimpleResidueMap& parser::residue::ResidueStorage::get(std::string name) {
    if (data.find(name) == data.end()) {
        std::cout << "Unknown residue: \"" << name << "\". Attempting to download specification." << std::endl;
        download_residue(name);
    }
    return data.at(name);
}

void parser::residue::ResidueStorage::initialize() {
    std::string path = setting::general::residue_folder;
    utility::create_directory(path);
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
            // std::cout << "Read residue " << residue << " from master file." << std::endl;

            // prepare map
            std::map<std::string, unsigned int> map;
            while (file.peek() != EOF) {
                std::getline(file, line);
                // stop if we reach the start of a new residue
                if (line.empty() || line.find("#") != std::string::npos) {
                    break;
                }

                // lines are of the form "atom hydrogens"
                std::vector<std::string> tokens = utility::split(line, " \n\r");
                std::string atom = tokens[0];
                unsigned int hydrogens = std::stoi(tokens[1]);
                map.emplace(atom, hydrogens);
            }
            insert(residue, std::move(map));
        }
    }
}

void parser::residue::ResidueStorage::download_residue(std::string name) {
    std::string path = setting::general::residue_folder;
    std::regex regex("[A-Z0-9]{2,3}");

    if (std::regex_match(name, regex)) {
        // check if the file already exists. if not, download it.
        if (!std::filesystem::exists(path + name + ".cif")) {
            curl::download("https://files.rcsb.org/ligands/view/" + name + ".cif", path + name + ".cif"); // download the cif file
        } else {
            std::cout << "Residue " << name << " is already downloaded, but not present in the master list. Reloading and adding it now." << std::endl;
        }

        // parse the cif file & add it to storage
        insert(name, parser::residue::detail::Residue::parse(path + name + ".cif").to_map());

        // write the residue to the master file
        write_residue(name);
    } else {
        throw except::map_error("Error in ResidueStorage::download_residue: Invalid residue name: \"" + name + "\"");
    }
}

void parser::residue::ResidueStorage::write_residue(std::string name) {
    std::string path = setting::general::residue_folder;
    utility::create_directory(path);
    std::ofstream file(path + "master.dat", std::ios::app); // open in append mode
    if (!file.is_open()) {throw except::io_error("Error in ResidueStorage::write_residue: Could not open file: " + path + "master" + ".dat");}

    // write the map to the master file
    auto map = get(name);
    file << "#" << "\n" << name << "\n";
    for (const auto& pair : map) {
        file << pair.first << " " << pair.second << "\n";
    }
    file << std::endl;
}

//### STREAM OPERATORS ###//
std::ostream& parser::residue::detail::operator<<(std::ostream& os, const parser::residue::detail::Atom& a) {os << a.to_string(); return os;}
std::ostream& parser::residue::detail::operator<<(std::ostream& os, const parser::residue::detail::Bond& b) {os << b.to_string(); return os;}
std::ostream& parser::residue::detail::operator<<(std::ostream& os, const parser::residue::detail::Residue& l) {os << l.to_string(); return os;}