#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <memory>

#include <utility/Constants.h>
#include <utility/Utility.h>
#include <utility/Exceptions.h>

#include <curl/curl.h>

template<typename K, typename V>
struct Storage {
    Storage() {}

    Storage(std::map<K, V> map) : data(map) {}

    V get(K key) {
        if (data.find(key) == data.end()) {
            throw except::map_error("Key " + key + " not found in map");
        }
        return data.at(key);
    }

    void insert(K key, V val) {
        data.emplace(key, val);
    }

    std::map<K, V> data;
};

struct Atom {
    Atom(std::string name, std::string symbol) : name(name), symbol(symbol) {
        valency = constants::valence::atomic.get(symbol);
    }

    void add_bond(std::string symbol, unsigned int order) {
        if (symbol == "H") {
            hydrogen_bonds++;
        }

        if (valency < order) {
            throw std::runtime_error("Atom " + name + " has no room for a bond of order " + std::to_string(order));
        }

        valency -= order;
    }

    std::string to_string() const {
        return "Atom " + name + " with valency " + std::to_string(valency) + " and " + std::to_string(hydrogen_bonds) + " hydrogen bonds";
    }

    friend std::ostream& operator<<(std::ostream& os, const Atom& a) {os << a.to_string(); return os;}

    std::string name, symbol;
    unsigned int valency;
    unsigned int hydrogen_bonds = 0;
};

struct Bond {
    Bond(std::string name1, std::string name2, unsigned int order) : name1(name1), name2(name2), order(order) {}

    std::string to_string() const {
        return "Bond " + name1 + (order == 1 ? " - " : " = ") + name2;
    }

    friend std::ostream& operator<<(std::ostream& os, const Bond& b) {os << b.to_string(); return os;}

    static unsigned int parse_order(std::string order) {
        if (order == "SING") {return 1;}
        else if (order == "DOUB") {return 2;}
        else {throw std::runtime_error("Invalid bond order: " + order);}
    }

    std::string name1, name2;
    unsigned int order;
};

class Ligand {
    public: 
        Ligand(std::string name) : name(name) {}

        Ligand(std::string name, std::vector<Atom> atoms, std::vector<Bond> bonds) : name(name), atoms(atoms) {
            apply_bond(bonds);
        }

        void add_atom(std::string name, std::string symbol) {
            name_map.insert({name, atoms.size()});
            atoms.push_back(Atom(name, symbol));
        }

        void apply_bond(const std::vector<Bond>& bonds) {
            for (const Bond& b : bonds) {
                apply_bond(b);
            }
        }

        void apply_bond(const Bond& bond) {
            Atom& a1 = atoms.at(name_map.at(bond.name1));
            Atom& a2 = atoms.at(name_map.at(bond.name2));
            a1.add_bond(a2.symbol, bond.order);
            a2.add_bond(a1.symbol, bond.order);
        }

        std::string to_string() const {
            std::stringstream ss;        
            for (const Atom& a : atoms) {
                ss << a << std::endl;
            }
            return ss.str();
        }

        Storage<std::string, unsigned int> to_map() const {
            Storage<std::string, unsigned int> map;
            for (const Atom& a : atoms) {
                map.insert(a.name, a.hydrogen_bonds);
            }
            return map;
        }

        friend std::ostream& operator<<(std::ostream& os, const Ligand& l) {os << l.to_string(); return os;}

    private: 
        std::string name;
        std::map<std::string, unsigned int> name_map;
        std::vector<Atom> atoms;        
};

Ligand parse(std::string filename) {
    std::ifstream file(filename);

    std::string line;
    Ligand ligand(utility::stem(filename));
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

#include <utility/Curl.h>
int main(int, char const**) {
    std::string file = "ASH";



    Ligand GLY = parse("GLY.cif");
    Ligand ASH = parse("ASH.cif");

    Storage glycine = std::map<std::string, unsigned int>{{"N", 1}, {"CA", 2}, {"C", 0}, {"O", 0}, {"OXT", 1}};
    Storage hydrogen_bonds = Storage<std::string, Storage<std::string, unsigned int>>{};

    hydrogen_bonds.insert("GLY", GLY.to_map());
    hydrogen_bonds.insert("ASH", ASH.to_map());    

    Storage<std::string, unsigned int> valence = std::map<std::string, unsigned int>{
        {"H", 1}, {"C", 4}, {"N", 3}, {"O", 2}, {"S", 2}, {"P", 3}, {"F", 1}, {"Cl", 1}, {"Br", 1}, {"I", 1}
    };

    std::cout << valence.get("H") << std::endl;
    std::cout << GLY << std::endl;

    return 0;
}