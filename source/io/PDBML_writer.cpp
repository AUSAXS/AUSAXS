#pragma once

// includes
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// my own includes
#include "Writer.h"

class PDBML_writer : public Writer {
public: 
    /** Constructor for the PDBML_writer class. 
     * @param filename the name of the input file
     */
    PDBML_writer(string filename) : Writer(filename) {};

    void write(vector<Atom*>* protein_atoms, vector<Atom*>* hydration_atoms) override {
        for (Atom* a : *protein_atoms) {file << to_pdbml(a) << endl;}
        for (Atom* a : *hydration_atoms) {file << to_pdbml(a) << endl;}
        return;
    }

private:
    /** Returns a standard multi-line PDBML format representation of this atom. 
     * @return a multi-line PDBML string representation of this atom. 
     */
    string to_pdbml(Atom* atom) {
        return (format("<PDBx:atom_site id=\"%1%\"> \
        \n    <PDBx:Cartn_x>%2%</PDBx:Cartn_x> \
        \n    <PDBx:Cartn_y>%3%</PDBx:Cartn_y> \
        \n    <PDBx:Cartn_z>%4%</PDBx:Cartn_z> \
        \n    <PDBx:occupancy>%5%</PDBx:occupancy> \
        \n    <PDBx:type_symbol>%6%</PDBx:type_symbol> \
        \n    <PDBx:label_comp_id>%7%</PDBx:label_comp_id> \
        \n</PDBx:atom_site>") % atom->get_serial() % atom->get_x() % atom->get_y() % atom->get_z() % atom->get_occupancy() % atom->get_symbol() % atom->get_comp()).str();
    }
};