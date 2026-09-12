#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/BodySplitter.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("BodySplitter::split by indices", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    SECTION("splits into correct number of bodies") {
        auto molecule = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
        REQUIRE(molecule.size_body() == 3);
    }

    SECTION("single split") {
        auto molecule = BodySplitter::split("tests/files/LAR1-2.pdb", {99});
        REQUIRE(molecule.size_body() == 2);
        CHECK(molecule.get_body(0).size_atom() > 0);
        CHECK(molecule.get_body(1).size_atom() > 0);
    }

    SECTION("atoms are preserved") {
        Molecule original("tests/files/LAR1-2.pdb");
        auto split = BodySplitter::split("tests/files/LAR1-2.pdb", {99});
        int total_atoms = 0;
        for (int i = 0; i < split.size_body(); ++i) {
            total_atoms += split.get_body(i).size_atom();
        }
        CHECK(total_atoms == original.get_body(0).size_atom());
    }
}

TEST_CASE("BodySplitter::split by chain", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    SECTION("splits by chain ID") {
        auto molecule = BodySplitter::split("tests/files/SASDJG5.pdb");
        REQUIRE(molecule.size_body() > 1);
        for (int i = 0; i < molecule.size_body(); ++i) {
            CHECK(molecule.get_body(i).size_atom() > 0);
        }
    }
}
