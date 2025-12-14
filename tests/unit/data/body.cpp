#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <form_factor/FormFactorType.h>
#include <constants/Constants.h>
#include <settings/MoleculeSettings.h>
#include <settings/GeneralSettings.h>
#include <math/Vector3.h>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("Body::Body") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;

    SECTION("default") {
        Body body;
        CHECK(body.size_atom() == 0);
        CHECK(body.size_water() == 0);
    }

    SECTION("vector<AtomFF>&") {
        std::vector<AtomFF> atoms = {
            AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
            AtomFF({4, 5, 6}, form_factor::form_factor_t::O)
        };
        Body body(atoms);
        REQUIRE(body.size_atom() == 2);
        CHECK(body.get_atom(0) == atoms[0]);
        CHECK(body.get_atom(1) == atoms[1]);
    }

    SECTION("vector<AtomFF>&&") {
        std::vector<AtomFF> atoms = {
            AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
            AtomFF({4, 5, 6}, form_factor::form_factor_t::O)
        };
        Body body(std::move(atoms));
        REQUIRE(body.size_atom() == 2);
        CHECK(body.get_atom(0).coordinates() == Vector3<double>{1, 2, 3});
        CHECK(body.get_atom(1).coordinates() == Vector3<double>{4, 5, 6});
    }

    SECTION("vector<Atom>&") {
        std::vector<Atom> atoms = {
            Atom({1, 2, 3}, 4),
            Atom({5, 6, 7}, 8)
        };
        Body body(atoms);
        REQUIRE(body.size_atom() == 2);
        CHECK(body.get_atom(0).coordinates() == atoms[0].coordinates());
        CHECK(body.get_atom(0).weight() == atoms[0].weight());
    }

    SECTION("vector<AtomFF>&, vector<Water>&") {
        std::vector<AtomFF> atoms = {
            AtomFF({1, 2, 3}, form_factor::form_factor_t::C)
        };
        std::vector<Water> waters = {
            Water({7, 8, 9})
        };
        Body body(atoms, waters);
        REQUIRE(body.size_atom() == 1);
        REQUIRE(body.size_water() == 1);
        CHECK(body.get_atom(0) == atoms[0]);
        CHECK(body.get_waters()->get()[0] == waters[0]);
    }

    SECTION("Body&") {
        std::vector<AtomFF> atoms = {
            AtomFF({1, 2, 3}, form_factor::form_factor_t::C)
        };
        Body body1(atoms);
        Body body2(body1);
        REQUIRE(body2.size_atom() == 1);
        CHECK(body2.get_atom(0) == atoms[0]);
        
        body2.get_atom(0).coordinates() = Vector3<double>{10, 20, 30};
        CHECK(body1.get_atom(0).coordinates() == Vector3<double>{1, 2, 3});
    }

    SECTION("Body&&") {
        std::vector<AtomFF> atoms = {
            AtomFF({1, 2, 3}, form_factor::form_factor_t::C)
        };
        Body body1(atoms);
        Body body2(std::move(body1));
        REQUIRE(body2.size_atom() == 1);
        CHECK(body2.get_atom(0) == atoms[0]);
    }
}

TEST_CASE("Body::get_atoms") {
    std::vector<AtomFF> atoms = {
        AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
        AtomFF({4, 5, 6}, form_factor::form_factor_t::O),
        AtomFF({7, 8, 9}, form_factor::form_factor_t::N)
    };
    Body body(atoms);

    SECTION("const") {
        const Body& cbody = body;
        const auto& result = cbody.get_atoms();
        REQUIRE(result.size() == 3);
        CHECK(result[0] == atoms[0]);
        CHECK(result[1] == atoms[1]);
        CHECK(result[2] == atoms[2]);
    }

    SECTION("non-const") {
        auto& result = body.get_atoms();
        REQUIRE(result.size() == 3);
        result[0].coordinates() = Vector3<double>{10, 20, 30};
        CHECK(body.get_atom(0).coordinates() == Vector3<double>{10, 20, 30});
    }
}

TEST_CASE("Body::get_atom") {
    std::vector<AtomFF> atoms = {
        AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
        AtomFF({4, 5, 6}, form_factor::form_factor_t::O)
    };
    Body body(atoms);

    SECTION("const") {
        const Body& cbody = body;
        CHECK(cbody.get_atom(0) == atoms[0]);
        CHECK(cbody.get_atom(1) == atoms[1]);
    }

    SECTION("non-const") {
        body.get_atom(0).coordinates() = Vector3<double>{10, 20, 30};
        CHECK(body.get_atom(0).coordinates() == Vector3<double>{10, 20, 30});
    }
}

TEST_CASE("Body::get_waters") {
    std::vector<Water> waters = {
        Water({1, 2, 3}),
        Water({4, 5, 6})
    };
    Body body(std::vector<AtomFF>{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, waters);

    SECTION("const") {
        const Body& cbody = body;
        const auto result = cbody.get_waters();
        REQUIRE(result->get().size() == 2);
        CHECK(result->get()[0] == waters[0]);
        CHECK(result->get()[1] == waters[1]);
    }

    SECTION("non-const") {
        auto result = body.get_waters();
        REQUIRE(result->get().size() == 2);
        result->get()[0].coordinates() = Vector3<double>{10, 20, 30};
        CHECK(body.get_waters()->get()[0].coordinates() == Vector3<double>{10, 20, 30});
    }
}

TEST_CASE("Body::get_cm") {
    SECTION("centered at origin") {
        std::vector<AtomFF> atoms = {
            AtomFF({-1, -1, -1}, form_factor::form_factor_t::C),
            AtomFF({-1,  1, -1}, form_factor::form_factor_t::C),
            AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C),
            AtomFF({ 1,  1, -1}, form_factor::form_factor_t::C),
            AtomFF({-1, -1,  1}, form_factor::form_factor_t::C),
            AtomFF({-1,  1,  1}, form_factor::form_factor_t::C),
            AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C),
            AtomFF({ 1,  1,  1}, form_factor::form_factor_t::C)
        };
        Body body(atoms);
        Vector3<double> cm = body.get_cm();
        CHECK(cm == Vector3<double>{0, 0, 0});
    }

    SECTION("offset") {
        std::vector<AtomFF> atoms = {
            AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
            AtomFF({3, 4, 5}, form_factor::form_factor_t::C)
        };
        Body body(atoms);
        Vector3<double> cm = body.get_cm();
        CHECK(cm == Vector3<double>{2, 3, 4});
    }
}

TEST_CASE("Body::get_molar_mass") {
    std::vector<AtomFF> atoms = {
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        AtomFF({1, 1, 1}, form_factor::form_factor_t::O),
        AtomFF({2, 2, 2}, form_factor::form_factor_t::N)
    };
    Body body(atoms);
    
    double expected = constants::mass::get_mass(form_factor::form_factor_t::C) * constants::Avogadro +
                      constants::mass::get_mass(form_factor::form_factor_t::O) * constants::Avogadro +
                      constants::mass::get_mass(form_factor::form_factor_t::N) * constants::Avogadro;
    CHECK_THAT(body.get_molar_mass(), Catch::Matchers::WithinRel(expected, 1e-6));
}

TEST_CASE("Body::get_absolute_mass") {
    std::vector<AtomFF> atoms = {
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        AtomFF({1, 1, 1}, form_factor::form_factor_t::O)
    };
    Body body(atoms);
    
    double expected = constants::mass::get_mass(form_factor::form_factor_t::C) +
                      constants::mass::get_mass(form_factor::form_factor_t::O);
    CHECK_THAT(body.get_absolute_mass(), Catch::Matchers::WithinRel(expected, 1e-6));
}

TEST_CASE("Body::get_total_atomic_charge") {
    std::vector<AtomFF> atoms = {
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        AtomFF({1, 1, 1}, form_factor::form_factor_t::O),
        AtomFF({2, 2, 2}, form_factor::form_factor_t::H)
    };
    Body body(atoms);
    
    double expected = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::C).I0() +
                      form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::O).I0() +
                      form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::H).I0();
    CHECK(body.get_total_atomic_charge() == expected);
}

TEST_CASE("Body::translate") {
    std::vector<AtomFF> atoms = {
        AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
        AtomFF({4, 5, 6}, form_factor::form_factor_t::O)
    };
    Body body(atoms);
    
    body.translate(Vector3<double>{10, 20, 30});
    CHECK(body.get_atom(0).coordinates() == Vector3<double>{11, 22, 33});
    CHECK(body.get_atom(1).coordinates() == Vector3<double>{14, 25, 36});
}

TEST_CASE("Body::size_atom") {
    SECTION("empty") {
        Body body;
        CHECK(body.size_atom() == 0);
    }

    SECTION("with atoms") {
        std::vector<AtomFF> atoms = {
            AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
            AtomFF({4, 5, 6}, form_factor::form_factor_t::O),
            AtomFF({7, 8, 9}, form_factor::form_factor_t::N)
        };
        Body body(atoms);
        CHECK(body.size_atom() == 3);
    }
}

TEST_CASE("Body::size_water") {
    SECTION("empty") {
        Body body;
        CHECK(body.size_water() == 0);
    }

    SECTION("with waters") {
        std::vector<Water> waters = {
            Water({1, 2, 3}),
            Water({4, 5, 6})
        };
        Body body(std::vector<AtomFF>{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, waters);
        CHECK(body.size_water() == 2);
    }
}

TEST_CASE("Body::clear_hydration") {
    std::vector<Water> waters = {
        Water({1, 2, 3}),
        Water({4, 5, 6})
    };
    Body body(std::vector<AtomFF>{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, waters);
    REQUIRE(body.size_water() == 2);
    
    body.clear_hydration();
    CHECK(body.size_water() == 0);
}

TEST_CASE("Body::equals_content") {
    std::vector<AtomFF> atoms = {
        AtomFF({1, 2, 3}, form_factor::form_factor_t::C)
    };
    
    SECTION("equal bodies") {
        Body body1(atoms);
        Body body2(atoms);
        CHECK(body1.equals_content(body2));
    }

    SECTION("different atoms") {
        Body body1(atoms);
        std::vector<AtomFF> atoms2 = {
            AtomFF({4, 5, 6}, form_factor::form_factor_t::O)
        };
        Body body2(atoms2);
        CHECK_FALSE(body1.equals_content(body2));
    }

    SECTION("different number of atoms") {
        Body body1(atoms);
        std::vector<AtomFF> atoms2 = {
            AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
            AtomFF({4, 5, 6}, form_factor::form_factor_t::O)
        };
        Body body2(atoms2);
        CHECK_FALSE(body1.equals_content(body2));
    }
}
