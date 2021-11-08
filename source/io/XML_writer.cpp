#pragma once

// includes
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <utility>
#include <filesystem>

// my own includes
#include "Writer.h"

using std::unique_ptr, std::shared_ptr;

class XML_writer : public Writer {
public: 
    /** Constructor for the XML_writer class. 
     * @param filename the name of the input file
     */
    XML_writer(string filename) : Writer(filename) {};

    void write(vector<shared_ptr<Atom>>* protein_atoms, vector<shared_ptr<Atom>>* hydration_atoms) override {
        std::filesystem::path p(filename);
        output << format("<PDBx:datablock datablockName=%1%>") % p.stem() << endl;
        output << "   <PDBx:atom_siteCategory>" << endl;
        for (auto const& a : *protein_atoms) {output << to_pdbml(a) << endl;}
        for (auto const& a : *hydration_atoms) {output << to_pdbml(a) << endl;}
        output << "   </PDBx:atom_siteCategory>" << endl;
        output << "</PDBx:datablock>" << endl;
        return;
    }

private:
    /** Returns a standard multi-line PDBML format representation of this atom. 
     * @return a multi-line PDBML string representation of this atom. 
     */
    string to_pdbml(shared_ptr<Atom> atom) {
        return (format("      <PDBx:atom_site id=\"%1%\"> \
        \n         <PDBx:Cartn_x>%2%</PDBx:Cartn_x> \
        \n         <PDBx:Cartn_y>%3%</PDBx:Cartn_y> \
        \n         <PDBx:Cartn_z>%4%</PDBx:Cartn_z> \
        \n         <PDBx:occupancy>%5%</PDBx:occupancy> \
        \n         <PDBx:type_symbol>%6%</PDBx:type_symbol> \
        \n         <PDBx:label_comp_id>%7%</PDBx:label_comp_id> \
        \n      </PDBx:atom_site>") % atom->get_serial() % atom->get_x() % atom->get_y() % atom->get_z() % atom->get_occupancy() % atom->get_symbol() % atom->get_comp()).str();
    }
};