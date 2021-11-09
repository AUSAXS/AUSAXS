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

const bool check_pdb = true;
const bool check_xml = false;

bool compareFiles(const std::string& p1, const std::string& p2) {
    std::ifstream f1(p1, std::ifstream::binary);
    std::ifstream f2(p2, std::ifstream::binary); 
    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }   

    string l1, l2;    
    Atom a1, a2;
    int max_i = 9999; 
    for (int i = 0; i < max_i; i++) {
        getline(f1, l1);
        getline(f2, l2);
        if (l1.empty() && l2.empty()) { // if both lines are empty, we're most likely at the end of both files
            return true;
        }

        // since a value of 5.90 is converted to 5.9 in the new file, we must manually compare entries where this can happen
        if (Record::get_type(l1.substr(0, 6)) == Record::ATOM) { 
            a1.parse_pdb(l1);
            a2.parse_pdb(l2);
            if (!(a1 == a2)) {
                print_err("File comparison failed on lines");
                cout << l1 << "|\n" << l2 << endl;
                return false;
            }
        } 

        // otherwise we just compare the lines themselves
        else {
            if (l1 != l2) {
                print_err("File comparison failed on lines");
                cout << l1 << "|\n" << l2 << endl;
                return false;
            }
        }
    }
    return false;
}

int main(int argc, char const *argv[])
{
    cout << "Testing io functionalities...\t\r" << std::flush;
    // test basic reading functionalities
    std::ofstream pdb_file("temp.pdb");
    std::ofstream xml_file("temp.xml");
    pdb_file << "ATOM      1  LEU HOH A 129         2.1     3.2     4.3  0.50 42.04           O  " << endl;
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
    if (check_xml) {
        Protein* protein = new Protein("temp.xml");
        protein->save("temp2.xml");
        protein = new Protein("temp2.xml");
        vector<shared_ptr<Atom>>* atoms = protein->get_hydration_atoms();
        shared_ptr<Atom> a = (*atoms)[0];

        // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
        // we now compare the loaded values with the expected.
        try {
            assert(a->get_serial() == 1);
            assert(a->get_x() == 2.1);
            assert(a->get_y() == 3.2);
            assert(a->get_z() == 4.3);
            assert(a->get_occupancy() == 0.50);
            assert(a->get_element() == "O");
            assert(a->get_resName() == "HOH");
        } catch (const std::exception& e) {
            print_err("XML input/output failed.");
        }
        remove("temp.xml");
        remove("temp2.xml");
    }

    if (check_pdb) {
        // check PDB io
        Protein* protein = new Protein("temp.pdb");
        protein->save("temp2.pdb");
        protein = new Protein("temp2.pdb");
        vector<shared_ptr<Atom>>* atoms = protein->get_hydration_atoms();
        shared_ptr<Atom> a = (*atoms)[0];

        // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
        // we now compare the loaded values with the expected.
        try {
            assert(a->get_serial() == 1);
            assert(a->get_x() == 2.1);
            assert(a->get_y() == 3.2);
            assert(a->get_z() == 4.3);
            assert(a->get_occupancy() == 0.50);
            assert(a->get_element() == "O");
            assert(a->get_resName() == "HOH");
        } catch (const std::exception& e) {
            print_err("PDB input/output failed.");
        }
        remove("temp.pdb");
        remove("temp2.pdb");
    }

    // load and copy each file in the data/ folder, and then compare the two files line-by-line
    // this is probably one of the strongest tests we can make for i/o
    for (const auto& file : std::filesystem::recursive_directory_iterator("data")) // loop over all files in the data/ directory
        if (check_pdb && file.path().extension() == ".pdb") { // check if the extension is .pdb
            Protein* protein = new Protein(file.path().string());
            protein->save("temp.pdb");
            try {
                assert(compareFiles(file.path().string(), "temp.pdb"));
            } catch (const std::exception& e) {
                print_err("PDB input/output failed for file " + file.path().string() + ".");
            }
            remove("temp.pdb");
        }

    cout << "\033[1;32m" << "All io tests passed." << "\033[0m" << endl;
    return 0;
}