#include "constants/ConstantsFwd.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <settings/All.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace data::record;

struct fixture {
    Atom a1  =  Atom(1, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1, -1, -1), 1, 0, constants::atom_t::C, "0");
    Water w1 = Water(1, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1, -1, -1), 1, 0, constants::atom_t::C, "0");
};

TEST_CASE_METHOD(fixture, "Water::Water") {
    SECTION("Atom&&") {
        Water w1(std::move(a1));
        CHECK(w1.serial == 1);
        CHECK(w1.name == "C");
        CHECK(w1.altLoc == "");
        CHECK(w1.resName == "HOH");
        CHECK(w1.chainID == 'A');
        CHECK(w1.resSeq == 1);
        CHECK(w1.iCode == "");
        CHECK(w1.coords == Vector3({-1, -1, -1}));
        CHECK(w1.occupancy == 1);
        CHECK(w1.tempFactor == 0);
        CHECK(w1.element == constants::atom_t::C);
        CHECK(w1.charge == "0");
        CHECK(w1.is_water() == true);
    }

    SECTION("Atom&") {
        Water w1(a1);
        CHECK(w1.serial == 1);
        CHECK(w1.name == "C");
        CHECK(w1.altLoc == "");
        CHECK(w1.resName == "HOH");
        CHECK(w1.chainID == 'A');
        CHECK(w1.resSeq == 1);
        CHECK(w1.iCode == "");
        CHECK(w1.coords == Vector3({-1, -1, -1}));
        CHECK(w1.occupancy == 1);
        CHECK(w1.tempFactor == 0);
        CHECK(w1.element == constants::atom_t::C);
        CHECK(w1.charge == "0");
        CHECK(w1.is_water() == true);
    }

    SECTION("Water&") {
        Water w2(w1);
        CHECK(w2.serial == 1);
        CHECK(w2.name == "C");
        CHECK(w2.altLoc == "");
        CHECK(w2.resName == "LYS");
        CHECK(w2.chainID == 'A');
        CHECK(w2.resSeq == 1);
        CHECK(w2.iCode == "");
        CHECK(w2.coords == Vector3({-1, -1, -1}));
        CHECK(w2.occupancy == 1);
        CHECK(w2.tempFactor == 0);
        CHECK(w2.element == constants::atom_t::C);
        CHECK(w2.charge == "0");
        CHECK(w2.is_water() == true);
    }
}

TEST_CASE_METHOD(fixture, "Water::get_mass") {
    CHECK(w1.get_mass() == constants::mass::get_mass(constants::atom_t::O) + 2*constants::mass::get_mass(constants::atom_t::H));
}

TEST_CASE_METHOD(fixture, "Water::get_type") {
    CHECK(w1.get_type() == RecordType::WATER);
}

TEST_CASE_METHOD(fixture, "Water::get_recName") {
    CHECK(w1.get_recName() == "HETATM");
}

TEST_CASE_METHOD(fixture, "Water::is_water") {
    CHECK(w1.is_water() == true);
}

TEST_CASE("Water::create_new_water") {
    Water w1 = Water::create_new_water(Vector3<double>({1, 2, 3}));
    CHECK(w1.serial == -1);
    CHECK(w1.name == "O");
    CHECK(w1.altLoc == "");
    CHECK(w1.resName == "HOH");
    CHECK(w1.chainID == ' ');
    CHECK(w1.resSeq == -1);
    CHECK(w1.iCode == "");
    CHECK(w1.coords == Vector3({1, 2, 3}));
    CHECK(w1.occupancy == 1);
    CHECK(w1.tempFactor == 0);
    CHECK(w1.element == constants::atom_t::O);
    CHECK(w1.charge == "");
    CHECK(w1.get_type() == RecordType::WATER);
    CHECK(w1.is_water() == true);
}

TEST_CASE("Water::operator=") {
    Water w2 = Water::create_new_water(Vector3<double>({1, 2, 3}));
    Water w3 = Water::create_new_water(Vector3<double>({1, 2, 3}));
    w2 = w3;
    CHECK(w2.serial == -1);
    CHECK(w2.name == "O");
    CHECK(w2.altLoc == "");
    CHECK(w2.resName == "HOH");
    CHECK(w2.chainID == ' ');
    CHECK(w2.resSeq == -1);
    CHECK(w2.iCode == "");
    CHECK(w2.coords == Vector3({1, 2, 3}));
    CHECK(w2.occupancy == 1);
    CHECK(w2.tempFactor == 0);
    CHECK(w2.element == constants::atom_t::O);
    CHECK(w2.charge == "");
    CHECK(w2.get_type() == RecordType::WATER);
    CHECK(w2.is_water() == true);
}

TEST_CASE("Water::operator==") {
    Water w2 = Water::create_new_water(Vector3<double>({1, 2, 3}));
    Water w3 = Water::create_new_water(Vector3<double>({1, 2, 3}));
    CHECK_FALSE(w2 == w3);

    w2 = w3;
    CHECK(w2 == w3);
}