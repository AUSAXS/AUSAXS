#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <fstream>
#include <iostream>
#include <memory>
#include <functional>

#include <utility/Constants.h>
#include <utility/Exceptions.h>

using std::string, std::cout, std::endl, std::vector;

class Ligand {
    class Atom {
        public: 
            /**
             * @brief Default constructor.
             */
            Atom() noexcept {}

            Atom(string symbol) {
                this->symbol = symbol;
                valence = constants::valence::get.at(symbol);
            }

            void add_connection(Ligand::Atom* atom) {
                if (valence <= 0 || atom->valence <= 0) {
                    throw std::runtime_error("Cannot add connections to an atom with no valence.");
                }
                connections.push_back(atom);
                atom->connections.push_back(this);
                atom->valence--;
                valence--;
            }
        
            string symbol, name;
            unsigned int valence;
            vector<Ligand::Atom*> connections;
    };

    public: 
        /**
         * @brief Default constructor.
         */
        Ligand() noexcept {}

        Ligand(string inchi) : inchi(inchi) {
            parse(inchi);
        }

        /**
         * @brief Parse the InChI string and create a ligand from it.
         */
        void parse(string id) {
            if (id.size() < 9 || id.substr(0, 9) != "InChI=1S/") {
                cout << "Fatal error, wrong prefix" << endl;
                exit(1);
            }

            unsigned int index = 9;
            std::map<unsigned int, Ligand::Atom> atom_ids = parse_initial(id, index);
            structure = parse_connections(id, index, atom_ids);
        }

        Ligand::Atom* find_calpha() const {
            for (unsigned int i = 0; i < structure.size(); i++) {
                Ligand::Atom atom1 = structure[i];

                // first search for the carboxylic acid group. we check all oxygens to see if they match the structure
                std::cout << atom1.symbol << std::endl;
                if (atom1.symbol == "O") {
                    std::cout << "Found oxygen" << std::endl;

                    // check if the oxygen is connected to a carbon
                    for (unsigned int j = 0; j < atom1.connections.size(); j++) {
                        Ligand::Atom* atom2 = atom1.connections[j];
                        if (atom2->symbol == "C") {
                            std::cout << "Found attached carbon" << std::endl;

                            // check if the carbon is connected to another oxygen
                            unsigned int oxygens = 0;
                            Ligand::Atom* carbon;
                            for (unsigned int k = 0; k < atom2->connections.size(); k++) {
                                std::cout << "k = " << k << std::endl;
                                Ligand::Atom* atom3 = atom2->connections[k];

                                // count the number of oxygens
                                if (atom3->symbol == "O") {
                                    std::cout << "Found attached oxygen" << std::endl;
                                    oxygens++;
                                } 

                                // save the connected carbon
                                else if (atom3->symbol == "C") {
                                    std::cout << "Found calpha" << std::endl;
                                    carbon = atom3;
                                }
                            }

                            // if there are two oxygens, we have found the carboxylic acid group and can easily find the calpha
                            if (oxygens == 2) {
                                carbon->name = "CA";
                                atom2->name = "C";
                                atom1.name = "OXT";
                                return carbon;
                            }
                        }
                    }
                }
            }
            return nullptr;
        }

        void enumerate() {
            Ligand::Atom* calpha = find_calpha();
            if (calpha == nullptr) {throw std::runtime_error("Cannot enumerate a ligand without a calpha.");}

            // explore a single atom and recursively explore its connections
            std::function<void(Ligand::Atom*, unsigned int)> explore = [&explore] (Ligand::Atom* atom, unsigned int level) {
                static std::map<unsigned int, char> level_key = {{0, 'A'}, {1, 'B'}, {2, 'G'}, {3, 'D'}, {4, 'E'}, {5, 'Z'}};

                atom->name = atom->symbol + std::to_string(level_key.at(level));
                for (unsigned int i = 0; i < atom->connections.size(); i++) {
                    if (atom->connections[i]->name.empty()) { // prevents backtracking and loops
                        explore(atom->connections[i], level + 1);
                    }
                }
            };

            // enumerate the atoms
            for (unsigned int i = 0; i < calpha->connections.size(); i++) {
                Ligand::Atom* atom = calpha->connections[i];
                if (atom->name == "C") {continue;}
                else if (atom->symbol == "N") {continue;}
                else {
                    explore(atom, 1);
                }
            }
        }

        void print() const noexcept {
            for (unsigned int i = 0; i < structure.size(); i++) {
                Ligand::Atom atom = structure[i];
                cout << atom.name << endl;
            }
        }

    private:
        string inchi;
        std::vector<Ligand::Atom> structure;

        /**
         * @brief Parse the chemical formula section of the InChI string.
         */
        std::map<unsigned int, Ligand::Atom> parse_initial(std::string id, unsigned int& index) {
            std::map<unsigned int, Ligand::Atom> atoms;
            unsigned int counter = 1;

            // iterate through the entire first section of the id
            while(index < id.size()) {
                //### determine the symbol ###//

                // check if we are at the end of the first section
                string symbol = string(1, id[index++]); // get the first character of the symbol
                if (symbol == "/") {
                    cout << "Found end of first section." << endl;
                    break;
                }

                // symbols can be either 1 or 2 characters long; the second character will always be lowercase
                if (std::islower(id[index])) {
                    symbol += id[index++];
                }
                cout << "symbol is " << symbol << endl;

                //### determine the number following the symbol ###//
                string number;
                // the number is optional; if it is not present, it is implicitly 1
                if (std::isdigit(id[index])) {
                    // while not at the end of the string
                    while(index < id.size()) {
                        // get the next character
                        char c = id[index];
                        // if it not a digit, we are at the end of the number
                        if (!std::isdigit(c)) {
                            break;
                        }
                        // add the character to the number
                        number += c;
                        index++;
                    }
                } else {
                    number = "1";
                }
                cout << "number is " << number << endl;
                // check we found a number
                if (number.empty()) {
                    cout << "Fatal error - invalid string format." << endl;
                    exit(2);
                }
                // check we're not dealing with the H (should be ignored here)
                if (symbol == "H") {
                    cout << "Found H, skipping" << endl;
                    continue;
                } 

                // everything is ok so we add to the map
                for (unsigned int i = 0; i < std::stoi(number); i++) {
                    atoms[counter++] = Ligand::Atom(symbol);
                }
            }

            return atoms;
        }

        /**
         * @brief Parse the connections section of the InChI string.
         */
        std::vector<Ligand::Atom> parse_connections(string id, unsigned int& index, std::map<unsigned int, Ligand::Atom>& atom_ids) {
            std::vector<Ligand::Atom> structure;

            if (id[index] != 'c') {
                std::cout << id[index] << std::endl;
                throw except::invalid_argument("Invalid InChI string format.");
            }
            index ++;

            while(index < id.size()) {
                // check if we are at the end of the second section
                if (id[index] == '/') {
                    cout << "Found end of second section." << endl;
                    break;
                }

                unsigned int current = 0, previous = 0;
                // parse the connection character between two atoms
                static auto parse_sep = [&previous, &current] (char sep) {
                    static std::vector<unsigned int> branchers;
                    switch (sep) {
                        case '-':
                            return true;
                        case '(':
                            branchers.push_back(current);
                            return true;
                        case ')':
                            if (branchers.size() == 0) {
                                throw except::invalid_argument("Invalid InChI string format.");
                            }
                            current = branchers.back();
                            branchers.pop_back();
                            return true;
                        case '/':
                            return false;
                        default: 
                            std::cout << "Invalid character: " << sep << std::endl;
                            throw except::invalid_argument("Invalid InChI string format.");
                    }
                };

                while(index < id.size()) {
                    string number;
                    while(index < id.size()) {
                        char c = id[index];
                        if (!std::isdigit(c)) {
                            break;
                        }
                        number += c;
                        index++;
                    }
                    previous = current;
                    current = std::stoi(number); 

                    // in the very first iteration no connection is made
                    if (previous == 0) {
                        parse_sep(id[index++]);
                        continue;
                    } 
                    // afterwards we always make at least one connection
                    else {
                        atom_ids[previous].add_connection(&atom_ids[current]);
                        if (!parse_sep(id[index++])) {
                            index--; // go back one character to avoid skipping the '/'
                            break;
                        }
                    }
                }
            }

            for (auto& atom : atom_ids) {structure.push_back(atom.second);}
            return structure;
        }
};

// int main(int, char const**) {
//     string test = "InChI=1S/C21H19ClN6O/c22-19-13-24-21-20(26-18(14-28(19)21)15-2-1-7-23-12-15)25-16-3-5-17(6-4-16)27-8-10-29-11-9-27/h1-7,12-14H,8-11H2,(H,25,26)";
//     Ligand ligand(test);
//     ligand.enumerate();
//     ligand.print();
//     return 0;
// }

#include <utility/Dataset.h>

int main(int, char const**) {
    return 0;
}