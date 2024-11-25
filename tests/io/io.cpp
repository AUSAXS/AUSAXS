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

using namespace ausaxs;
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
                std::cout << l1 << "|\n" << l2 << "|" << std::endl;
                return false;
            }
        }

        // sometimes nothing is written after TER in the pdb files
        else if (type1 == RecordType::TERMINATE) {continue;}

        // otherwise we just compare the lines themselves
        else {
            if (l1 != l2) {
                console::print_warning("File line comparison failed for \"" + p1 + "\" on lines");
                std::cout << l1 << "|\n" << l2 << "|" << std::endl;
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
    pdb_file << "REMARK ONE" << std::endl;
    pdb_file << "REMARK TWO" << std::endl;
    pdb_file << "CRYST1 THREE" << std::endl;
    pdb_file << "ATOM      1  CB  ARG A 129         2.1     3.2     4.3  0.50 42.04           C " << std::endl;
    pdb_file << "TER       2      ARG A 129                                                     " << std::endl;
    pdb_file << "HETATM    3  O   HOH A 130      31.117   3.049  35.879  0.94 34.19           O " << std::endl;
    pdb_file << "HETATM    4  O   HOH A 131      31.117   3.049  35.879  0.94 34.19           O " << std::endl;
    pdb_file << "MASTER FOUR" << std::endl;
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

TEST_CASE("io: protein") {
    settings::molecule::center = false;
    settings::general::verbose = false;

    Molecule protein("tests/files/2epe.pdb");
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
    Body body("tests/files/2epe.pdb");
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