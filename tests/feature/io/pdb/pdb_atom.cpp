#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <cmath>

#include <io/pdb/PDBAtom.h>
#include <io/pdb/PDBWater.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

TEST_CASE("PDBAtom::distance") {
    SECTION("simple") {
        PDBAtom a1({0, 0, 0}, 1, constants::atom_t::N, "GLY", 1);
        PDBAtom a2({1, 0, 0}, 2, constants::atom_t::C, "GLY", 1);
        CHECK(a1.coordinates().distance(a2.coordinates()) == 1);
    }

    SECTION("complex") {
        PDBAtom a1({0, 0, 0}, 1, constants::atom_t::N, "GLY", 1);
        PDBAtom a2({1, 0, 0}, 2, constants::atom_t::C, "GLY", 1);
        PDBAtom a3({1, 1, 0}, 3, constants::atom_t::O, "GLY", 1);
        CHECK(a1.coordinates().distance(a2.coordinates()) == 1);
        CHECK(a1.coordinates().distance(a3.coordinates()) == std::sqrt(2));
    }
}

TEST_CASE("PDBAtom::translate") {
    SECTION("simple") {
        PDBAtom a1({0, 0, 0}, 1, constants::atom_t::N, "GLY", 1);
        a1.coordinates() += Vector3{1, 2, 3};
        CHECK(a1.coordinates() == Vector3{1, 2, 3});
    }

    SECTION("complex") {
        PDBAtom a1({0, 0, 0}, 1, constants::atom_t::N, "GLY", 1);
        a1.coordinates() += Vector3{1, 2, 3};
        CHECK(a1.coordinates() == Vector3<double>{1, 2, 3});
        a1.coordinates() += Vector3{1, 2, 3};
        CHECK(a1.coordinates() == Vector3<double>{2, 4, 6});
    }
}

TEST_CASE("PDBAtom: operators") {
    PDBAtom a1({3, 0, 5}, 2, constants::atom_t::He, "", 3);
    PDBAtom a2 = a1;
    REQUIRE(a1 == a2);
    REQUIRE(a1.equals_content(a2));

    a2 = PDBAtom({0, 4, 1}, 2, constants::atom_t::He, "", 2);
    REQUIRE(a1 != a2);
    REQUIRE(!a1.equals_content(a2));
    REQUIRE(a2 < a1);

    PDBWater w1 = PDBWater({3, 0, 5}, 2, constants::atom_t::He, "", 3);
    PDBWater w2 = w1;
    REQUIRE(w1 == w2);

    w2 = PDBAtom({0, 4, 1}, 2, constants::atom_t::He, "", 2);
    REQUIRE(w1 != w2);
    REQUIRE(w2 < w1);
}
