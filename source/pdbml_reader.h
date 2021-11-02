#pragma once

// includes
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// my own includes
#include "Reader.h"

class pdbml_reader : Reader {
public: 

    /** Constructor for the pdbml_reader class. 
     * @param filename the name of the input file
     */
    pdbml_reader(std::string filename) : Reader(filename) {
        if (filename.find(".xml") != std::string::npos) {
            perror(("Input file \"" + filename + "\" is not a .xml file!").c_str());
            exit(EXIT_FAILURE);
        }
    };

    std::vector<Atom> read() override {
        std::string line;
        Atom* atom = new Atom();
        std::vector<Atom*> atoms;

        while(getline(file, line)) {
            // check if this line contains any relevant information
            if (line.find("<PDBx:") == std::string::npos) {
                continue; // otherwise we just skip it
            }

            // check if the line is the start of an atom_site
            if (line.find("<PDBx:atom_site>") != std::string::npos) {
                // since this is the start of an atom_site, we make a sanity check on the serial number (it should be the previous serial+1)
                int v_start = line.find("\"") + 1;
                int v_end = line.find_last_of("\"");
                std::string v = line.substr(v_start, v_end);

                // sanity check
                if (atoi(v.c_str()) != atom->serial+1) {
                    perror(("Broken reading sequence after " + v + ". Terminating.").c_str());
                    exit(EXIT_FAILURE);
                } else {
                    atom->serial = atoi(v.c_str());
                }
                continue;
            }

            // check if the line is the end of an atom site
            else if (line.find("</PDBx:atom_site>") != std::string::npos) {
                // at this point, the atom should contain all of the relevant properties, so we simply store it in an array and prepare another one.
                atoms.push_back(atom);
                atom = new Atom();
                continue;
            }

            // the remaining options all follow the same format of "<PDBx:(w)>(v)</PDBx:(w)>". we want to find (w) and (v)
            int w_start = line.find("<PDBx:"); // the characters in front of w
            int w_end = line.find_first_of(">"); // the character after w
            int v_start = w_end+1; // the start of v
            int v_end = line.find("</PDBx:"); // the characters after v

            std::string w = line.substr(w_start, w_end);
            std::string v = line.substr(v_start, v_end);

            // debug
            std::cout << "word is: " << w << std::endl;
            std::cout << "value is: " << v << std::endl;
            
            // check if the word is a coordinate
            if (w == "Cartn_x") {
                atom->x = stod(v);
                continue;
            } else if (w == "Cartn_y") {
                atom->y = stod(v);
                continue;
            } else if (w == "Cartn_z") {
                atom->z = stod(v);
                continue;
            }
            
            // check if the word is the weight
            else if (w == "occupancy") {
                atom->w = stod(v);
                continue;
            }
        }
    };

private:
    Atom read_line() override {

    };
};