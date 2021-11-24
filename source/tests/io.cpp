// includes
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

#include "tests/Test.h"
#include "Protein.h"

using namespace ROOT;

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
        Record::RecordType type = Record::get_type(l1.substr(0, 6));
        if (type == Record::ATOM || type == Record::HETATM) { 
            a1.parse_pdb(l1);
            a2.parse_pdb(l2);
            if (!(a1 == a2)) {
                print_err("File comparison failed for \"" + p1 + "\" on lines");
                cout << l1 << "|\n" << l2 << "|" << endl;
                return false;
            }
        } 

        // otherwise we just compare the lines themselves
        else {
            if (l1 != l2) {
                print_err("File comparison failed for \"" + p1 + "\" on lines");
                cout << l1 << "|\n" << l2 << "|" << endl;
                return false;
            }
        }
    }
    return false;
}

/**
 * @brief Perform a simple unit test of .pdb files. 
 */
void test_simple_pdb() {
    std::ofstream pdb_file("temp.pdb");
    pdb_file << "ATOM      1  LEU   O A 129         2.1     3.2     4.3  0.50 42.04           O " << endl;
    pdb_file << "HETATM    2  C1  MYR A 601      31.117   3.049  35.879  0.94 34.19           C " << endl;
    pdb_file << "HETATM    3  HOH HOH A 601      31.117   3.049  35.879  0.94 34.19           O " << endl;
    pdb_file.close();

    // check PDB io
    Protein* protein = new Protein("temp.pdb");
    protein->save("temp2.pdb");
    protein = new Protein("temp2.pdb");
    vector<shared_ptr<Atom>>* atoms = protein->get_protein_atoms();
    shared_ptr<Atom> a = (*atoms)[0];

    if (atoms->size() == 0) {
        IS_TRUE(false);
        return;
    }

    // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
    // we now compare the loaded values with the expected.
    IS_TRUE(a->get_serial() == 1);
    IS_TRUE(a->get_x() == 2.1);
    IS_TRUE(a->get_y() == 3.2);
    IS_TRUE(a->get_z() == 4.3);
    IS_TRUE(a->get_occupancy() == 0.50);
    IS_TRUE(a->get_element() == "O");
    IS_TRUE(a->get_resName() == "O");

    remove("temp.pdb");
    remove("temp2.pdb");
}

/**
 * @brief Perform a simple io test of .xml files.
 */
void test_simple_xml() {
    std::ofstream xml_file("temp.xml");
    xml_file << "<PDBx:atom_site id=\"1\"> \
        \n    <PDBx:Cartn_x>2.1</PDBx:Cartn_x> \
        \n    <PDBx:Cartn_y>3.2</PDBx:Cartn_y> \
        \n    <PDBx:Cartn_z>4.3</PDBx:Cartn_z> \
        \n    <PDBx:occupancy>0.50</PDBx:occupancy> \
        \n    <PDBx:type_symbol>O</PDBx:type_symbol> \
        \n    <PDBx:label_comp_id>HOH</PDBx:label_comp_id> \
        \n</PDBx:atom_site>" << endl;
    xml_file.close();

    Protein* protein = new Protein("temp.xml");
    protein->save("temp2.xml");
    protein = new Protein("temp2.xml");
    vector<shared_ptr<Hetatom>>* atoms = protein->get_hydration_atoms();
    shared_ptr<Atom> a = (*atoms)[0];

    // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
    // we now compare the loaded values with the expected.
    IS_TRUE(a->get_serial() == 1) 
    IS_TRUE(a->get_serial() == 1);
    IS_TRUE(a->get_x() == 2.1);
    IS_TRUE(a->get_y() == 3.2);
    IS_TRUE(a->get_z() == 4.3);
    IS_TRUE(a->get_occupancy() == 0.50);
    IS_TRUE(a->get_element() == "O");
    IS_TRUE(a->get_resName() == "HOH");

    remove("temp.xml");
    remove("temp2.xml");
}

/**
 * @brief Load and copy each file in the data/ folder, and then compare the two files line-by-line.
 *        This is probably one of the strongest tests we can make for i/o
 */
void test_all_data() {
    for (const auto& file : std::filesystem::recursive_directory_iterator("data")) { // loop over all files in the data/ directory
        if (file.path().extension() == ".pdb") { // check if the extension is .pdb
            Protein* protein = new Protein(file.path().string());
            protein->save("temp.pdb");
            IS_TRUE(compareFiles(file.path().string(), "temp.pdb"));
            remove("temp.pdb");
        }
    }
}

int main(void) {
    cout << "Summary of IO testing:" << endl;
    print_err("Test that TER records are generated correctly.");

    test_simple_pdb();
    // test_simple_xml();
    test_all_data();

    if (passed_all) {
        cout << "\033[1;32m" << "All io tests passed.           " << "\033[0m" << endl;
    } else {
        print_err("Some IO tests failed.");
    }
    return 0;
}