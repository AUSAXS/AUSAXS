#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <data/atoms/Atom.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <form_factor/FormFactorType.h>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("Atom::Atom") {
    SECTION("default") {
        Atom atom;
        CHECK_THAT(atom.weight(), Catch::Matchers::WithinAbs(0, 1e-6));
    }

    SECTION("Vector3<precision_t>&, precision_t") {
        Vector3 coords = GENERATE(
            Vector3<double>{1, 2, 3},
            Vector3<double>{4, 5, 6},
            Vector3<double>{7, 8, 9}
        );
        double w = GENERATE(1., 2., 3.);

        Atom atom(coords, w);
        CHECK(atom.coordinates() == coords);
        CHECK(atom.weight() == w);
        CHECK(atom.x() == coords.x());
        CHECK(atom.y() == coords.y());
        CHECK(atom.z() == coords.z());
    }
}

TEST_CASE("Atom::coordinates") {
    Atom atom({1, 2, 3}, 4);
    
    SECTION("const getter") {
        CHECK(atom.coordinates() == Vector3<double>{1, 2, 3});
        CHECK(atom.position() == Vector3<double>{1, 2, 3});
    }

    SECTION("setter") {
        atom.coordinates() = Vector3<double>{5, 6, 7};
        CHECK(atom.coordinates() == Vector3<double>{5, 6, 7});
    }

    SECTION("individual accessors") {
        CHECK(atom.x() == 1);
        CHECK(atom.y() == 2);
        CHECK(atom.z() == 3);
        
        atom.x() = 10;
        atom.y() = 20;
        atom.z() = 30;
        CHECK(atom.coordinates() == Vector3<double>{10, 20, 30});
    }
}

TEST_CASE("Atom::weight") {
    Atom atom({1, 2, 3}, 4);
    
    SECTION("const getter") {
        CHECK(atom.weight() == 4);
    }

    SECTION("setter") {
        atom.weight() = 5;
        CHECK(atom.weight() == 5);
    }
}

TEST_CASE("Atom::equality") {
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

TEST_CASE("AtomFF::AtomFF") {
    SECTION("default") {
        AtomFF atom;
        CHECK_THAT(atom.weight(), Catch::Matchers::WithinAbs(0, 1e-6));
    }

    SECTION("Atom&, form_factor_t") {
        Atom basic({1, 2, 3}, 4);
        AtomFF atom(basic, form_factor::form_factor_t::C);
        CHECK(atom.coordinates() == basic.coordinates());
        CHECK(atom.weight() == basic.weight());
        CHECK(atom.form_factor_type() == form_factor::form_factor_t::C);
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
        AtomFF atom({1, 2, 3}, form_factor::form_factor_t::C, 0.5);
        CHECK(atom.coordinates() == Vector3<double>{1, 2, 3});
        CHECK(atom.form_factor_type() == form_factor::form_factor_t::C);
        CHECK(atom.weight() == 0.5);
    }
}

TEST_CASE("AtomFF::form_factor_type") {
    AtomFF atom({1, 2, 3}, form_factor::form_factor_t::C);
    
    SECTION("const getter") {
        CHECK(atom.form_factor_type() == form_factor::form_factor_t::C);
    }

    SECTION("setter") {
        atom.form_factor_type() = form_factor::form_factor_t::O;
        CHECK(atom.form_factor_type() == form_factor::form_factor_t::O);
    }
}

TEST_CASE("AtomFF::equality") {
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

TEST_CASE("Water::Water") {
    SECTION("default") {
        Water water;
        CHECK_THAT(water.weight(), Catch::Matchers::WithinAbs(0, 1e-6));
    }

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

TEST_CASE("Water::coordinates") {
    Water water({1, 2, 3});
    
    SECTION("const getter") {
        CHECK(water.coordinates() == Vector3<double>{1, 2, 3});
        CHECK(water.position() == Vector3<double>{1, 2, 3});
    }

    SECTION("setter") {
        water.coordinates() = Vector3<double>{5, 6, 7};
        CHECK(water.coordinates() == Vector3<double>{5, 6, 7});
    }
}

TEST_CASE("Water::equality") {
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
