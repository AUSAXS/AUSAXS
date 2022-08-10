#include <utility/ResidueParser.h>
#include <utility/Curl.h>
#include <utility/Settings.h>
#include <utility/Constants.h>

#include <regex>
#include <filesystem>

//### ATOM ###//
parser::ligand::detail::Atom::Atom(std::string name, std::string symbol) : name(name), symbol(symbol) {
    valency = constants::valence::atomic.get(symbol);
}

void parser::ligand::detail::Atom::add_bond(std::string symbol, unsigned int order) {
    if (symbol == "H") {
        hydrogen_bonds++;
    }

    if (valency < order) {
        throw std::runtime_error("Atom " + name + " has no room for a bond of order " + std::to_string(order));
    }

    valency -= order;
}

std::string parser::ligand::detail::Atom::to_string() const {
    return "Atom " + name + " with valency " + std::to_string(valency) + " and " + std::to_string(hydrogen_bonds) + " hydrogen bonds";
}


//### BOND ###//
parser::ligand::detail::Bond::Bond(std::string name1, std::string name2, unsigned int order) : name1(name1), name2(name2), order(order) {}

std::string parser::ligand::detail::Bond::to_string() const {
    return "Bond " + name1 + (order == 1 ? " - " : " = ") + name2;
}

unsigned int parser::ligand::detail::Bond::parse_order(std::string order) {
    if (order == "SING") {return 1;}
    else if (order == "DOUB") {return 2;}
    else {throw std::runtime_error("Invalid bond order: " + order);}
}


//### RESIDUE ###//
parser::ligand::detail::Residue::Residue(std::string name) : name(name) {}

parser::ligand::detail::Residue::Residue(std::string name, std::vector<Atom> atoms, std::vector<Bond> bonds) : name(name), atoms(atoms) {
    apply_bond(bonds);
}

void parser::ligand::detail::Residue::add_atom(std::string name, std::string symbol) {
    name_map.insert({name, atoms.size()});
    atoms.push_back(Atom(name, symbol));
}

void parser::ligand::detail::Residue::apply_bond(const std::vector<Bond>& bonds) {
    for (const Bond& b : bonds) {
        apply_bond(b);
    }
}

void parser::ligand::detail::Residue::apply_bond(const Bond& bond) {
    Atom& a1 = atoms.at(name_map.at(bond.name1));
    Atom& a2 = atoms.at(name_map.at(bond.name2));
    a1.add_bond(a2.symbol, bond.order);
    a2.add_bond(a1.symbol, bond.order);
}

std::string parser::ligand::detail::Residue::to_string() const {
    std::stringstream ss;        
    for (const Atom& a : atoms) {
        ss << a << std::endl;
    }
    return ss.str();
}

std::map<std::string, unsigned int> parser::ligand::detail::Residue::to_map() const {
    std::map<std::string, unsigned int> map;
    for (const Atom& a : atoms) {
        map.emplace(a.name, a.hydrogen_bonds);
    }
    return map;
}

parser::ligand::detail::Residue parser::ligand::detail::Residue::parse(std::string filename) {
    std::ifstream file(filename);

    std::string line;
    Residue ligand(utility::stem(filename));
    bool found_atom_section = false, found_bond_section = false;
    while (std::getline(file, line)) {
        if (line.find("atom.pdbx_ordinal") != std::string::npos) {
            std::cout << "Found start of atom section" << std::endl;
            found_atom_section = true;
            break;
        }
    }

    while (std::getline(file, line)) {
        // check for end of section
        if (line.find("#") != std::string::npos) {
            std::cout << "Found end of atom section" << std::endl;
            break;
        }

        std::vector<std::string> tokens = utility::split(line, " \n\r");
        std::string atom_id = tokens[1];
        std::string type_symbol = tokens[3];
        ligand.add_atom(atom_id, type_symbol);
    }

    while (std::getline(file, line)) {
        if (line.find("bond.pdbx_ordinal") != std::string::npos) {
            std::cout << "Found start of bond section" << std::endl;
            found_bond_section = true;
            break;
        }
    }

    while (std::getline(file, line)) {
        // check for end of section
        if (line.find("#") != std::string::npos) {
            std::cout << "Found end of bond section" << std::endl;
            break;
        }

        std::vector<std::string> tokens = utility::split(line, " \n\r");
        std::string atom1 = tokens[1];
        std::string atom2 = tokens[2];
        unsigned int order = Bond::parse_order(tokens[3]);
        ligand.apply_bond(Bond(atom1, atom2, order));
    }
    
    return ligand;
}

//### RESIDUE STORAGE ###//
parser::ligand::ResidueStorage::ResidueStorage() {
    initialize();
}

void parser::ligand::ResidueStorage::insert(std::string name, std::map<std::string, unsigned int> ligand) {
    data.emplace(name, ligand);
}

std::map<std::string, unsigned int>& parser::ligand::ResidueStorage::get(std::string name) {
    if (data.find(name) == data.end()) {
        download_residue(name);

    }
    return data.at(name);
}

void parser::ligand::ResidueStorage::initialize() {
    std::string path = setting::general::residue_folder;
    // utility::create_directories(path);

    for (const auto& file : std::filesystem::directory_iterator(path)) {
        if (file.path().extension() == ".cif") {
            std::string name = file.path().string();
            insert(name, parser::ligand::detail::Residue::parse(path + name + ".cif").to_map());
        }
    }
}

void parser::ligand::ResidueStorage::download_residue(std::string name) {
    std::string path = setting::general::residue_folder;
    std::regex regex("[A-Z]{3}");

    if (std::regex_match(name, regex)) {
        curl::download("https://files.rcsb.org/ligands/view/" + name + ".cif", path + name + ".cif");
        insert(name, parser::ligand::detail::Residue::parse(path + name + ".cif").to_map());
    } else {
        throw except::map_error("Error in ResidueStorage::download_ligand: Invalid ligand name: \"" + name + "\"");
    }
}

//### STREAM OPERATORS ###//
std::ostream& parser::ligand::detail::operator<<(std::ostream& os, const parser::ligand::detail::Atom& a) {os << a.to_string(); return os;}
std::ostream& parser::ligand::detail::operator<<(std::ostream& os, const parser::ligand::detail::Bond& b) {os << b.to_string(); return os;}
std::ostream& parser::ligand::detail::operator<<(std::ostream& os, const parser::ligand::detail::Residue& l) {os << l.to_string(); return os;}