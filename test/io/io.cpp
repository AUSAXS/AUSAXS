#include "settings/MoleculeSettings.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Water.h>
#include <data/record/Record.h>
#include <utility/Console.h>
#include <settings/All.h>

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>

using std::cout, std::endl, std::vector;

using namespace data;
using namespace data::record;

bool compare_files(std::string p1, std::string p2) {
    std::ifstream f1(p1, std::ifstream::binary);
    std::ifstream f2(p2, std::ifstream::binary); 
    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }   

    std::string l1, l2;
    Atom a1, a2;
    int max_i = 99999; 
    for (int i = 0; i < max_i; i++) {
        getline(f1, l1);
        getline(f2, l2);
        if (l1.empty()) {
            if (l2.empty()) {return true;} // if both lines are empty, we're at the end of both files
            if (Record::get_type(l2.substr(0, 6)) == RecordType::TERMINATE) {return true;} // we allow a single terminate of difference
            console::print_warning("File ended prematurely.");
            return false;
        }

        RecordType type1 = Record::get_type(l1.substr(0, 6)); 
        RecordType type2 = Record::get_type(l2.substr(0, 6)); 
        if (type1 != type2) {
            console::print_warning("The types " + l1.substr(0, 6) + " and " + l2.substr(0, 6) + " are not equal in line " + std::to_string(i) + ".");
            return false;
        }

        // since a value of 5.90 is converted to 5.9 in the new file, we must manually compare entries where this can happen
        if (type1 == RecordType::ATOM) { 
            a1.parse_pdb(l1);
            a2.parse_pdb(l2);

            // equality of atoms is based on their unique ID which is generated at object creation. Thus this will never be equal with this approach.
            // instead we must compare their contents. 
            if (!a1.equals_content(a2)) {
                console::print_warning("File atom comparison failed for \"" + p1 + "\" on lines");
                cout << l1 << "|\n" << l2 << "|" << endl;
                return false;
            }
        }

        // sometimes nothing is written after TER in the pdb files
        else if (type1 == RecordType::TERMINATE) {continue;}

        // otherwise we just compare the lines themselves
        else {
            if (l1 != l2) {
                console::print_warning("File line comparison failed for \"" + p1 + "\" on lines");
                cout << l1 << "|\n" << l2 << "|" << endl;
                return false;
            }
        }
    }
    return false;
}

TEST_CASE("io: body file") {
    settings::general::verbose = false;
    io::File path("temp/io/temp.pdb");
    path.create();

    std::ofstream pdb_file(path);
    pdb_file << "REMARK ONE" << endl;
    pdb_file << "REMARK TWO" << endl;
    pdb_file << "CRYST1 THREE" << endl;
    pdb_file << "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C " << endl;
    pdb_file << "TER       2      ARG A 129                                                     " << endl;
    pdb_file << "HETATM    3  O   HOH A 130      31.117   3.049  35.879  0.94 34.19           O " << endl;
    pdb_file << "HETATM    4  O   HOH A 131      31.117   3.049  35.879  0.94 34.19           O " << endl;
    pdb_file << "MASTER FOUR" << endl;
    pdb_file.close();

    Body body(path);
    auto& file = body.get_file();
    REQUIRE(file.header.size() == 3);
    REQUIRE(file.footer.size() == 1);

    file.header.remove("CRYST1");
    REQUIRE(file.header.size() == 2);
    REQUIRE(file.header.get() == "REMARK ONE\nREMARK TWO\n");

    auto& file2 = body.get_file();
    REQUIRE(file2.header.size() == 2);

    path = "temp/io/temp2.pdb";
    body.save("temp/io/temp2.pdb");
    Body body2("temp/io/temp2.pdb");
    auto file3 = body2.get_file();
    REQUIRE(file3.header.size() == 2);
}

TEST_CASE("io: pdb input") {
    settings::molecule::center = false;
    settings::general::verbose = false;
    settings::grid::scaling = 2;

    io::File path("temp/io/temp.pdb");
    path.create();

    std::ofstream pdb_file(path);
    pdb_file << "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C " << endl;
    pdb_file << "ATOM      2  CB  ARG A 129         3.2     4.3     5.4  0.50 42.04           C " << endl;
    pdb_file << "TER       3      ARG A 129                                                     " << endl;
    pdb_file << "HETATM    4  O   HOH A 130      30.117  29.049  34.879  0.94 34.19           O " << endl;
    pdb_file << "HETATM    5  O   HOH A 131      31.117  30.049  35.879  0.94 34.19           O " << endl;
    pdb_file.close();

    // check PDB io
    std::unique_ptr<Molecule> protein = std::make_unique<Molecule>(path);
    path = "temp/io/temp2.pdb";
    protein->save(path);
    protein = std::make_unique<Molecule>(path);
    vector<Atom> atoms = protein->get_atoms();
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
    CHECK(a.element == constants::atom_t::C);
    CHECK(a.resName == "ARG");
}

TEST_CASE("io: xml input", "[broken]") {
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

    Molecule* protein = new Molecule("temp.xml");
    protein->save("temp2.xml");
    protein = new Molecule("temp2.xml");
    const vector<Water>& atoms = protein->get_waters();
    const Water a = atoms[0];

    // the idea is that we have now loaded the hardcoded strings above, saved them, and loaded them again. 
    // we now compare the loaded values with the expected.
    CHECK(a.serial == 1);
    CHECK(a.serial == 1);
    CHECK(a.coords.x() == 2.1);
    CHECK(a.coords.y() == 3.2);
    CHECK(a.coords.z() == 4.3);
    CHECK(a.occupancy == 0.50);
    CHECK(a.element == constants::atom_t::O);
    CHECK(a.resName == "HOH");
}

/**
 * @brief Load and copy each file in the data/ folder, and then compare the two files line-by-line.
 *        This is probably one of the strongest tests we can make for i/o
 */
TEST_CASE("io: real data", "[files],[broken]") {
    settings::general::verbose = false;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    for (const auto& file : std::filesystem::recursive_directory_iterator("data")) { // loop over all files in the data/ directory
        if (file.path().extension() != ".pdb") {
            continue;
        }

        if (file.path().string().find("urateox") != std::string::npos) { // skip this file since it skips some serials (shouldn't this be illegal?????)
            cout << "Skipped urateox.pdb" << endl;
            continue;
        }

        cout << "Testing " << file.path().stem() << endl;
        std::string filename = "temp/io/" + file.path().stem().string() + ".pdb";
        Molecule protein(file.path().string());
        protein.save(filename);
        bool success = compare_files(file.path().string(), filename);
        CHECK(success);
    }
}

TEST_CASE("io: protein") {
    settings::molecule::center = false;
    settings::molecule::use_effective_charge = false;
    settings::general::verbose = false;

    Molecule protein("test/files/2epe.pdb");
    protein.save("temp/io/temp.pdb");
    Molecule protein2("temp/io/temp.pdb");
    auto atoms1 = protein.get_atoms();
    auto atoms2 = protein2.get_atoms();

    REQUIRE(atoms1.size() == atoms2.size());
    for (unsigned int i = 0; i < atoms1.size(); i++) {
        REQUIRE(atoms1[i].equals_content(atoms2[i]));
    }

    auto waters1 = protein.get_waters();
    auto waters2 = protein2.get_waters();
    REQUIRE(waters1.size() == waters2.size());
    for (unsigned int i = 0; i < waters1.size(); i++) {
        waters2[i].set_chainID(waters1[i].get_chainID()); // we always use a new chainID for the hydration
        bool equal = waters1[i].equals_content(waters2[i]);
        if (!equal) {
            std::cout << waters1[i].as_pdb() << std::endl;
            std::cout << waters2[i].as_pdb() << std::endl;
        }
        REQUIRE(equal);
    }
}

TEST_CASE("io: body copying") {
    Body body("test/files/2epe.pdb");
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

TEST_CASE("io: writing multifile pdb") {
    settings::molecule::implicit_hydrogens = false;

    std::vector<Atom> atoms(101000);
    for (unsigned int i = 0; i < atoms.size(); i++) {
        atoms[i] = Atom({1,2,3}, 1, constants::atom_t::C, "LYS", i);
    }
    std::vector<Water> waters(100);
    for (unsigned int i = 0; i < waters.size(); i++) {
        waters[i] = Water({1,2,3}, 1, constants::atom_t::O, "HOH", i);
    }

    Molecule protein(atoms, waters);
    protein.save("temp/io/temp_multifile.pdb");
    
    REQUIRE(std::filesystem::exists("temp/io/temp_multifile_part1.pdb"));
    REQUIRE(std::filesystem::exists("temp/io/temp_multifile_part2.pdb"));

    // first file
    Molecule protein2("temp/io/temp_multifile_part1.pdb");
    REQUIRE(protein2.get_body(0).atom_size() == 100000);
    REQUIRE(protein2.get_body(0).get_atoms().back().serial == 99999);

    Molecule protein3("temp/io/temp_multifile_part2.pdb");
    REQUIRE(protein3.get_body(0).atom_size() == 1000);
    for (int i = 0; i < 1000; i++) {
        REQUIRE(protein3.get_body(0).get_atom(i).serial == i);
    }
    REQUIRE(protein3.get_body(0).get_file().terminate.serial == 1000);
    REQUIRE(protein3.water_size() == 100);
    for (int i = 0; i < 100; i++) {
        REQUIRE(protein3.get_water(i).serial == i+1001);
    }
}

#include <em/ImageStack.h>
TEST_CASE("io: writing em multifile pdb", "[broken]") {
    em::ImageStack image("data/Gregers_cryo/Gregers_cryo.mrc");
    image.get_protein(0.1)->save("temp/io/temp.pdb");
}

TEST_CASE("can_parse_hydrogens") {
    std::vector<std::string> val = {"ATOM      1  N   VAL     1      -3.299   8.066 -11.443  1.00  0.00           N",
                                    "ATOM      2  H   VAL     1      -3.411   8.677 -12.239  1.00  0.00           H",
                                    "ATOM      3  CA  VAL     1      -3.085   6.673 -11.780  1.00  0.00           C",
                                    "ATOM      4  HA  VAL     1      -3.328   6.080 -10.899  1.00  0.00           H",
                                    "ATOM      5  CB  VAL     1      -3.927   6.165 -12.947  1.00  0.00           C",
                                    "ATOM      6  HB  VAL     1      -3.774   6.930 -13.708  1.00  0.00           H",
                                    "ATOM      7  CG1 VAL     1      -3.577   4.780 -13.486  1.00  0.00           C",
                                    "ATOM      8 HG11 VAL     1      -3.508   4.011 -12.716  1.00  0.00           H",
                                    "ATOM      9 HG12 VAL     1      -4.289   4.438 -14.237  1.00  0.00           H",
                                    "ATOM     10 HG13 VAL     1      -2.612   4.760 -13.992  1.00  0.00           H",
                                    "ATOM     11  CG2 VAL     1      -5.370   6.065 -12.463  1.00  0.00           C",
                                    "ATOM     12 HG21 VAL     1      -5.876   7.021 -12.328  1.00  0.00           H",
                                    "ATOM     13 HG22 VAL     1      -6.043   5.567 -13.162  1.00  0.00           H",
                                    "ATOM     14 HG23 VAL     1      -5.354   5.510 -11.524  1.00  0.00           H",
                                    "ATOM     15  C   VAL     1      -1.621   6.436 -12.123  1.00  0.00           C",
                                    "ATOM     16  O   VAL     1      -1.200   7.045 -13.104  1.00  0.00           O",
    };

    Atom atom;
    atom.parse_pdb(val[7]);
    REQUIRE(atom.name == "HG11");
    REQUIRE(atom.element == constants::atom_t::H);
    REQUIRE_THAT(atom.coords.x(), Catch::Matchers::WithinAbs(-3.508, 1e-6));
    REQUIRE_THAT(atom.coords.y(), Catch::Matchers::WithinAbs(4.011, 1e-6));
    REQUIRE_THAT(atom.coords.z(), Catch::Matchers::WithinAbs(-12.716, 1e-6));

    atom.parse_pdb(val[8]);
    REQUIRE(atom.name == "HG12");
    REQUIRE(atom.element == constants::atom_t::H);
    REQUIRE_THAT(atom.coords.x(), Catch::Matchers::WithinAbs(-4.289, 1e-6));
    REQUIRE_THAT(atom.coords.y(), Catch::Matchers::WithinAbs(4.438, 1e-6));
    REQUIRE_THAT(atom.coords.z(), Catch::Matchers::WithinAbs(-14.237, 1e-6));

    atom.parse_pdb(val[9]);
    REQUIRE(atom.name == "HG13");
    REQUIRE(atom.element == constants::atom_t::H);
    REQUIRE_THAT(atom.coords.x(), Catch::Matchers::WithinAbs(-2.612, 1e-6));
    REQUIRE_THAT(atom.coords.y(), Catch::Matchers::WithinAbs(4.760, 1e-6));
    REQUIRE_THAT(atom.coords.z(), Catch::Matchers::WithinAbs(-13.992, 1e-6));
}