#pragma once

// includes
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// my own includes
#include "Reader.h"
#include "Atom.cpp"

class pdbml_reader : public Reader {
public: 

    /** Constructor for the pdbml_reader class. 
     * @param filename the name of the input file
     */
    pdbml_reader(std::string filename) : Reader(filename) {
        if (filename.find(".xml") == std::string::npos) {
            perror(("Input file \"" + filename + "\" is not a .xml file!").c_str());
            exit(EXIT_FAILURE);
        }
    };

    std::vector<Atom*> read() override {
        std::string line;
        Atom* atom = new Atom();
        std::vector<Atom*> atoms;
        int serial = 1;
        atom->set_serial(serial);

        while(getline(file, line)) {
            // check if this line contains any relevant information
            if (line.find("PDBx:") == std::string::npos) {
                continue; // otherwise we just skip it
            }

            // check if the line is the start of an atom_site
            if (line.find("<PDBx:atom_site ") != std::string::npos) {
                // since this is the start of an atom_site, we make a sanity check on the serial number (it should be the previous serial+1)
                int v_start = line.find("\"")+1;
                int v_end = line.find_last_of("\"");
                std::string v = line.substr(v_start, v_end-v_start);

                // sanity check
                if (atoi(v.c_str()) != serial) {
                    perror(("Broken reading sequence after " + std::to_string(serial) + ". Terminating.").c_str());
                    exit(EXIT_FAILURE);
                }
                continue;
            }

            // check if the line is the end of an atom site
            else if (line.find("</PDBx:atom_site>") != std::string::npos) {
                // at this point, the atom should contain all of the relevant properties, so we simply store it in an array and prepare another one.
                atoms.push_back(atom);
                atom = new Atom();

                serial++;
                atom->set_serial(serial);
                continue;
            }

            // the remaining options all follow the same format of "<PDBx:(w)>(v)</PDBx:(w)>". we want to find (w) and (v)
            int w_start = line.find("<PDBx:")+6; // the characters in front of w
            int w_end = line.find_first_of(">"); // the character after w
            int v_start = w_end+1; // the start of v
            int v_end = line.find("</PDBx:"); // the characters after v

            std::string w = line.substr(w_start, w_end-w_start);
            std::string v = line.substr(v_start, v_end-v_start);

            // debug
            // std::cout << "word is: (" << w << ")" << std::endl;
            // std::cout << "value is: (" << v << ")" << std::endl;
            
            // check if the word is a coordinate
            if (w == "Cartn_x") {
                atom->set_x(stod(v));
                continue;
            } else if (w == "Cartn_y") {
                atom->set_y(stod(v));
                continue;
            } else if (w == "Cartn_z") {
                atom->set_z(stod(v));
                continue;
            }
            
            // check if the word is the weight
            else if (w == "occupancy") {
                atom->set_weight(stod(v));
                continue;
            }

            // check if the word is the symbol
            else if (w == "type_symbol") {
                atom->set_symbol(v);
                continue;
            }

            // check if the word describes the molecule
            else if (w == "auth_comp_id") {
                atom->set_comp(v);
                continue;
            }

            std::cout << "no match found" << std::endl;
        }
        return atoms;
    };

private:
    Atom read_line() override {
        return Atom();
    };
};