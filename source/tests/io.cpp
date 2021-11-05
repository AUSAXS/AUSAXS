// includes
#include <vector>
#include <string>
#include <fstream>

// ROOT
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1.h>

// my own includes
#include "../Protein.cpp"

using namespace ROOT;

int main(int argc, char const *argv[])
{
    cout << "Testing io functionalities...\t\r" << std::flush;
    // test basic reading functionalities
    std::ofstream pdb_file("temp.pdb");
    std::ofstream xml_file("temp.xml");
    pdb_file << "ATOM      1  LEU HOH A 129       2.1     3.2     4.3    0.50 42.04           O  " << endl;
    xml_file << "<PDBx:atom_site id=\"1\"> \
        \n    <PDBx:Cartn_x>2.1</PDBx:Cartn_x> \
        \n    <PDBx:Cartn_y>3.2</PDBx:Cartn_y> \
        \n    <PDBx:Cartn_z>4.3</PDBx:Cartn_z> \
        \n    <PDBx:occupancy>0.50</PDBx:occupancy> \
        \n    <PDBx:type_symbol>O</PDBx:type_symbol> \
        \n    <PDBx:label_comp_id>HOH</PDBx:label_comp_id> \
        \n</PDBx:atom_site>" << endl;
    pdb_file.close();
    xml_file.close();

    // check XML io
    Protein* protein = new Protein("temp.xml");
    protein->save("temp2.xml");
    protein = new Protein("temp2.xml");
    vector<Atom*>* atoms = protein->get_hydration_atoms();
    Atom* a = (*atoms)[0];

    // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
    // we now compare the loaded values with the expected.
    try {
        assert(a->get_serial() == 1);
        assert(a->get_x() == 2.1);
        assert(a->get_y() == 3.2);
        assert(a->get_z() == 4.3);
        assert(a->get_occupancy() == 0.50);
        assert(a->get_symbol() == "O");
        assert(a->get_comp() == "HOH");
    } catch (const std::exception& e) {
        print_err("XML input/output failed.");
    }
    remove("temp.xml");
    remove("temp2.xml");

    // check PDB io
    protein = new Protein("temp.pdb");
    protein->save("temp2.pdb");
    protein = new Protein("temp2.pdb");
    atoms = protein->get_hydration_atoms();
    a = (*atoms)[0];

    // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
    // we now compare the loaded values with the expected.
    try {
        assert(a->get_serial() == 1);
        assert(a->get_x() == 2.1);
        assert(a->get_y() == 3.2);
        assert(a->get_z() == 4.3);
        assert(a->get_occupancy() == 0.50);
        assert(a->get_symbol() == "O");
        assert(a->get_comp() == "HOH");
    } catch (const std::exception& e) {
        print_err("PDB input/output failed.");
    }
    remove("temp.pdb");
    remove("temp2.pdb");

    cout << "\033[1;32m" << "All io tests passed." << "\033[0m" << endl;
    return 0;
}