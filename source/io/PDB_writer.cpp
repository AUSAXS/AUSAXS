#pragma once

// includes
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// my own includes
#include "Writer.h"

using std::setw, std::left, std::right;

class PDB_writer : public Writer {
public: 
    /** Constructor for the PDBML_writer class. 
     * @param filename the name of the input file
     */
    PDB_writer(string filename) : Writer(filename) {};

    void write(vector<Atom*>* protein_atoms, vector<Atom*>* hydration_atoms) override {
        auto write_v = [&] (vector<Atom*>* atoms) {
            for (Atom* a : *atoms) {
                file << left << setw(6) << "ATOM" << " "            // starts at index 0
                << right << setw(4) << a->get_serial() << "  "      // 7
                << left << setw(3) << "NAN" << " "                  // 13
                << left << setw(3) << a->get_comp() << " "          // 17
                << "A" << " "                                       // 21
                << right << setw(3) << "00" << " "                  // 23
                << "    "                                           // 27
                << right << setw(7) << a->get_x() << " "            // 31
                << right << setw(7) << a->get_y() << " "            // 39
                << right << setw(7) << a->get_z() << " "            // 47
                << right << setw(5) << a->get_occupancy() << " "    // 55
                << right << setw(5) << "00.00" << " "               // 61
                << "        "                                       // 67
                << right << setw(3) << a->get_symbol() << " "       // 75
                << endl;                                            // 79
            }
        };
        write_v(protein_atoms);
        file << left << setw(6) << "TER" << " " << right << setw(4) << protein_atoms->size() << "  " << endl;
        write_v(hydration_atoms);
        return;
    }
};