#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>

#include <data/Protein.h>
#include <utility/Utility.h>

using std::cout, std::endl;

bool compareFiles(const std::string& p1, const std::string& p2) {
    std::ifstream f1(p1, std::ifstream::binary);
    std::ifstream f2(p2, std::ifstream::binary); 
    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }   

    string l1, l2;
    Atom a1, a2;
    int max_i = 99999; 
    for (int i = 0; i < max_i; i++) {
        getline(f1, l1);
        getline(f2, l2);
        if (l1.empty()) {
            if (l2.empty()) {return true;} // if both lines are empty, we're at the end of both files
            if (Record::get_type(l2.substr(0, 6)) == Record::TERMINATE) {return true;} // we allow a single terminate of difference
            utility::print_warning("File ended prematurely.");
            return false;
        }

        Record::RecordType type1 = Record::get_type(l1.substr(0, 6)); 
        Record::RecordType type2 = Record::get_type(l2.substr(0, 6)); 
        if (type1 != type2) {
            utility::print_warning("The types " + l1.substr(0, 6) + " and " + l2.substr(0, 6) + " are not equal in line " + std::to_string(i) + ".");
            return false;
        }

        // since a value of 5.90 is converted to 5.9 in the new file, we must manually compare entries where this can happen
        if (type1 == Record::ATOM || type1 == Record::HETATM) { 
            a1.parse_pdb(l1);
            a2.parse_pdb(l2);

            // equality of atoms is based on their unique ID which is generated at object creation. Thus this will never be equal with this approach.
            // instead we must compare their contents. 
            if (!a1.equals_content(a2)) {
                utility::print_warning("File atom comparison failed for \"" + p1 + "\" on lines");
                cout << l1 << "|\n" << l2 << "|" << endl;
                return false;
            }
        }

        // sometimes nothing is written after TER in the pdb files
        else if (type1 == Record::TERMINATE) {continue;}

        // otherwise we just compare the lines themselves
        else {
            if (l1 != l2) {
                utility::print_warning("File line comparison failed for \"" + p1 + "\" on lines");
                cout << l1 << "|\n" << l2 << "|" << endl;
                return false;
            }
        }
    }
    return false;
}

TEST_CASE("body_file", "[io]") {
    std::ofstream pdb_file("temp/io/temp.pdb");
    pdb_file << "REMARK ONE" << endl;
    pdb_file << "REMARK TWO" << endl;
    pdb_file << "CRYST1 THREE" << endl;
    pdb_file << "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C " << endl;
    pdb_file << "TER       2      ARG A 129                                                     " << endl;
    pdb_file << "HETATM    3  O   HOH A 130      31.117   3.049  35.879  0.94 34.19           O " << endl;
    pdb_file << "HETATM    4  O   HOH A 131      31.117   3.049  35.879  0.94 34.19           O " << endl;
    pdb_file << "MASTER FOUR" << endl;
    pdb_file.close();

    Body body("temp/io/temp.pdb");
    auto& file = body.get_file();
    REQUIRE(file.header.size() == 3);
    REQUIRE(file.footer.size() == 1);

    file.header.remove("CRYST1");
    REQUIRE(file.header.size() == 2);
    REQUIRE(file.header.get() == "REMARK ONE\nREMARK TWO\n");

    auto& file2 = body.get_file();
    REQUIRE(file2.header.size() == 2);

    body.save("temp/io/temp2.pdb");
    Body body2("temp/io/temp2.pdb");
    auto file3 = body2.get_file();
    REQUIRE(file3.header.size() == 2);

    remove("temp/io/temp.pdb");
    remove("temp/io/temp2.pdb");
}

TEST_CASE("pdb_input", "[io]") {
    std::ofstream pdb_file("temp/io/temp.pdb");
    pdb_file << "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C " << endl;
    pdb_file << "ATOM      2  CB  ARG A 129         3.2     4.3     5.4  0.50 42.04           C " << endl;
    pdb_file << "TER       3      ARG A 129                                                     " << endl;
    pdb_file << "HETATM    4  O   HOH A 130      31.117   3.049  35.879  0.94 34.19           O " << endl;
    pdb_file << "HETATM    5  O   HOH A 131      31.117   3.049  35.879  0.94 34.19           O " << endl;
    pdb_file.close();

    // check PDB io
    Protein* protein = new Protein("temp/io/temp.pdb");
    protein->save("temp/io/temp2.pdb");
    protein = new Protein("temp/io/temp2.pdb");
    vector<Atom> atoms = protein->get_protein_atoms();
    Atom a = atoms[0];

    if (atoms.size() == 0) {
        REQUIRE(false);
        return;
    }

    // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
    // we now compare the loaded values with the expected.
    CHECK(a.serial == 1);
    CHECK(a.coords.x() == 2.1);
    CHECK(a.coords.y() == 3.2);
    CHECK(a.coords.z() == 4.3);
    CHECK(a.occupancy == 0.50);
    CHECK(a.element == "C");
    CHECK(a.resName == "ARG");

    remove("temp/io/temp.pdb");
    remove("temp/io/temp2.pdb");
}

TEST_CASE("xml input", "[broken],[io]") {
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
    const vector<Hetatom>& atoms = protein->get_hydration_atoms();
    const Atom a = atoms[0];

    // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
    // we now compare the loaded values with the expected.
    CHECK(a.serial == 1);
    CHECK(a.serial == 1);
    CHECK(a.coords.x() == 2.1);
    CHECK(a.coords.y() == 3.2);
    CHECK(a.coords.z() == 4.3);
    CHECK(a.occupancy == 0.50);
    CHECK(a.element == "O");
    CHECK(a.resName == "HOH");

    remove("temp.xml");
    remove("temp2.xml");
}

/**
 * @brief Load and copy each file in the data/ folder, and then compare the two files line-by-line.
 *        This is probably one of the strongest tests we can make for i/o
 */
TEST_CASE("real_data", "[io],[files]") {
    for (const auto& file : std::filesystem::recursive_directory_iterator("data")) { // loop over all files in the data/ directory
        if (file.path() == "data/6yg9.pdb") { // skip this file since it contains OQ5 ligands which we can't deal with yet
            cout << "Skipped 6yg9.pdb" << endl;
            continue;
        }
        if (file.path() == "data/urateox.pdb") { // skip this file since it skips some serials (shouldn't this be illegal?????)
            cout << "Skipped urateox.pdb" << endl;
            continue;
        }

        if (file.path().extension() == ".pdb") { // check if the extension is .pdb
            cout << "Testing " << file.path().stem() << endl;
            Protein protein(file.path().string());
            protein.save("temp.pdb");
            REQUIRE(compareFiles(file.path().string(), "temp.pdb"));
            remove("temp.pdb");
        }
    }
}

TEST_CASE("file_copied_correctly", "[io],[files]") {
    Body body("data/lysozyme/2epe.pdb");
    CHECK(!body.get_file().header.get().empty());
    CHECK(!body.get_file().footer.get().empty());

    Body body2 = body;
    CHECK(!body2.get_file().header.get().empty());
    CHECK(!body2.get_file().footer.get().empty());

    Body body3;
    body3 = body;
    CHECK(!body3.get_file().header.get().empty());
    CHECK(!body3.get_file().footer.get().empty());
}