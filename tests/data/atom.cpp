#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <data/atoms/Atom.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <form_factor/FormFactorType.h>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("Atom::Atom") {
    SECTION("Vector3<precision_t>&, precision_t") {
        Vector3 coords = GENERATE(
            Vector3<double>{1, 2, 3},
            Vector3<double>{4, 5, 6},
            Vector3<double>{7, 8, 9}
        );
        double w = GENERATE(1., 2., 3.);

        Atom atom(coords, w);
        CHECK(atom.coordinates() == coords);
        CHECK(atom.position() == coords);
        CHECK(atom.x() == coords.x());
        CHECK(atom.y() == coords.y());
        CHECK(atom.z() == coords.z());
        CHECK(atom.weight() == w);
    }
}

TEST_CASE("AtomFF::AtomFF") {
    SECTION("const Atom&, form_factor_t") {
        Vector3 coords = GENERATE(
            Vector3<double>{1, 2, 3},
            Vector3<double>{4, 5, 6},
            Vector3<double>{7, 8, 9}
        );
        form_factor::form_factor_t ff = GENERATE(
            form_factor::form_factor_t::C, 
            form_factor::form_factor_t::O, 
            form_factor::form_factor_t::H
        );
        Atom atom(coords, constants::charge::nuclear::get_charge(ff));
        AtomFF atom_ff(atom, ff);
        CHECK(atom_ff.coordinates() == coords);
        CHECK(atom_ff.form_factor_type() == ff);
        CHECK(atom_ff.weight() == constants::charge::nuclear::get_charge(ff));
        CHECK(atom_ff.get_atom() == atom);
    }

    SECTION("Vector3<precision_t>&, form_factor_t") {
        Vector3 coords = GENERATE(
            Vector3<double>{1, 2, 3},
            Vector3<double>{4, 5, 6},
            Vector3<double>{7, 8, 9}
        );
        form_factor::form_factor_t ff = GENERATE(
            form_factor::form_factor_t::C, 
            form_factor::form_factor_t::O, 
            form_factor::form_factor_t::H
        );
        AtomFF atom(coords, ff);
        CHECK(atom.coordinates() == coords);
        CHECK(atom.form_factor_type() == ff);
        CHECK(atom.weight() == constants::charge::nuclear::get_charge(ff));
    }

    SECTION("Vector3<precision_t>&, form_factor_t, double") {
        Vector3 coords = GENERATE(
            Vector3<double>{1, 2, 3},
            Vector3<double>{4, 5, 6},
            Vector3<double>{7, 8, 9}
        );
        form_factor::form_factor_t ff = GENERATE(
            form_factor::form_factor_t::C, 
            form_factor::form_factor_t::O, 
            form_factor::form_factor_t::H
        );
        double w = GENERATE(1., 2., 3.);
        AtomFF atom(coords, ff, w);
        CHECK(atom.coordinates() == coords);
        CHECK(atom.form_factor_type() == ff);
        CHECK(atom.weight() == w);
    }
}

TEST_CASE("Water::Water") {
    SECTION("Vector3<precision_t>&") {
        Vector3 coords = GENERATE(
            Vector3<double>{1, 2, 3},
            Vector3<double>{4, 5, 6},
            Vector3<double>{7, 8, 9}
        );
        Water water(coords);
        CHECK(water.coordinates() == coords);
        CHECK(water.weight() == constants::charge::nuclear::get_charge(form_factor::form_factor_t::OH));
        CHECK(water.form_factor_type() == form_factor::form_factor_t::OH);
    }
}

TEST_CASE("Atom::operator==") {
    Atom a1({1, 2, 3}, 4);
    Atom a2({1, 2, 3}, 4);
    Atom a3({1, 2, 3}, 5);
    Atom a4({1, 2, 4}, 4);
    Atom a5({1, 3, 3}, 4);
    Atom a6({2, 2, 3}, 4);

    CHECK(a1 == a2);
    CHECK_FALSE(a1 == a3);
    CHECK_FALSE(a1 == a4);
    CHECK_FALSE(a1 == a5);
    CHECK_FALSE(a1 == a6);
}

TEST_CASE("AtomFF::operator==") {
    AtomFF a1({1, 2, 3}, form_factor::form_factor_t::C);
    AtomFF a2({1, 2, 3}, form_factor::form_factor_t::C);
    AtomFF a3({1, 2, 3}, form_factor::form_factor_t::O);
    AtomFF a4({1, 2, 4}, form_factor::form_factor_t::C);
    AtomFF a5({1, 3, 3}, form_factor::form_factor_t::C);
    AtomFF a6({2, 2, 3}, form_factor::form_factor_t::C);

    CHECK(a1 == a2);
    CHECK_FALSE(a1 == a3);
    CHECK_FALSE(a1 == a4);
    CHECK_FALSE(a1 == a5);
    CHECK_FALSE(a1 == a6);
}

TEST_CASE("Water::operator==") {
    Water w1({1, 2, 3});
    Water w2({1, 2, 3});
    Water w3({1, 2, 4});
    Water w4({1, 3, 3});
    Water w5({2, 2, 3});

    CHECK(w1 == w2);
    CHECK_FALSE(w1 == w3);
    CHECK_FALSE(w1 == w4);
    CHECK_FALSE(w1 == w5);
}