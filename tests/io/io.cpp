#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <io/pdb/PDBStructure.h>
#include <utility/Console.h>
#include <settings/All.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::io::pdb;

bool compare_files(std::string p1, std::string p2) {
    std::ifstream f1(p1, std::ifstream::binary);
    std::ifstream f2(p2, std::ifstream::binary); 
    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }   

    std::string l1, l2;
    PDBAtom a1, a2;
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
        REQUIRE(atoms1[i] == atoms2[i]);
    }

    auto waters1 = protein.get_waters();
    auto waters2 = protein2.get_waters();
    REQUIRE(waters1.size() == waters2.size());
    for (unsigned int i = 0; i < waters1.size(); i++) {
        REQUIRE(waters1[i] == waters2[i]);
    }
}