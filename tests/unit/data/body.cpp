#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <data/state/BoundSignaller.h>
#include <data/state/StateManager.h>
#include <form_factor/FormFactorInfo.h>
#include <hist/histogram_manager/IPartialHistogramManager.h>
#include <io/ExistingFile.h>
#include <io/Writer.h>
#include <math/Vector3.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <numbers>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::data;

struct fixture {
    std::vector<AtomFF> a = {
        AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C),
        AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C),
        AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C),
        AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)
    };
    Body body = Body(a);
};

struct multiple_fixture {
    AtomFF a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3 = AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a4 = AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
    AtomFF a5 = AtomFF( {1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a6 = AtomFF( {1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a7 = AtomFF( {1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a8 = AtomFF( {1,  1,  1}, form_factor::form_factor_t::H);
    Water w1 = Water({0, 1, 2});
    Water w2 = Water({3, 4, 5});

    Body b1 = Body(std::vector<AtomFF>{a1, a2});
    Body b2 = Body(std::vector<AtomFF>{a3, a4});
    Body b3 = Body(std::vector<AtomFF>{a5, a6});
    Body b4 = Body(std::vector<AtomFF>{a7, a8});
    std::vector<Body> ap = {b1, b2, b3, b4};
};

TEST_CASE("Body::Body") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

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

    SECTION("ExistingFile&") {
        settings::general::verbose = false;
        io::ExistingFile file("tests/files/2epe.pdb");
        Body b(file);
        CHECK(b.size_atom() == 1001);
        CHECK(b.size_water() == 48);
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
    
    // 1 Da is 1 g/mol by definition, so the molar mass in g/mol is numerically the summed atomic mass in Da. 
    double expected = constants::mass::get_mass(form_factor::form_factor_t::C) +
                      constants::mass::get_mass(form_factor::form_factor_t::O) +
                      constants::mass::get_mass(form_factor::form_factor_t::N);
    CHECK_THAT(body.get_molar_mass(), Catch::Matchers::WithinRel(expected, 1e-5));

    // independent order-of-magnitude anchor: C + N + O is ~42 g/mol
    CHECK_THAT(body.get_molar_mass(), Catch::Matchers::WithinAbs(42.02, 0.5));
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
    
    double expected = constants::charge::get_ff_charge(form_factor::form_factor_t::C) +
                      constants::charge::get_ff_charge(form_factor::form_factor_t::O) +
                      constants::charge::get_ff_charge(form_factor::form_factor_t::H);
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

TEST_CASE("Body::save") {
    settings::molecule::implicit_hydrogens = false;
    settings::general::verbose = false;
    std::vector<AtomFF> a = {
        AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C),
        AtomFF({ 1, -1, -1}, form_factor::form_factor_t::O), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C),
        AtomFF({-1, -1,  1}, form_factor::form_factor_t::N), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C),
        AtomFF({ 1, -1,  1}, form_factor::form_factor_t::O), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)
    };
    Body body(a);
    test::TempFile path(".pdb");
    io::Writer::write({body}, path);
    Body body2(path);

    CHECK(body.size_atom() == body2.size_atom());
    for (int i = 0; i < body.size_atom(); i++) {
        CHECK(body.get_atom(i) == body2.get_atom(i));
    }
}

TEST_CASE_METHOD(fixture, "Body::translate: informs the histogram manager") {
    auto protein = Molecule({body});
    protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManager);
    auto* manager = static_cast<hist::IPartialHistogramManager*>(protein.get_histogram_manager())->get_state_manager();
    manager->reset_to_false();
    protein.get_body(0).translate(Vector3<double>(10, 0, 0));
    CHECK(protein.get_body(0).get_atom(0).coordinates() == Vector3<double>(9, -1, -1));
    CHECK(protein.get_body(0).get_atom(1).coordinates() == Vector3<double>(9,  1, -1));
    CHECK(manager->get_externally_modified_bodies()[0] == true);
}

TEST_CASE("Body::rotate") {
    SECTION("simple") {
        std::vector<AtomFF> a = {
            AtomFF({1, 0, 0}, form_factor::form_factor_t::C),
            AtomFF({0, 1, 0}, form_factor::form_factor_t::C),
            AtomFF({0, 0, 1}, form_factor::form_factor_t::C)};
        Body body(a);

        Vector3<double> axis = {0, 1, 0};
        body.rotate(matrix::rotation_matrix(axis, std::numbers::pi/2));
        CHECK(Vector3<double>{0, 0, -1} == body.get_atom(0).coordinates());
        CHECK(Vector3<double>{0, 1,  0} == body.get_atom(1).coordinates());
        CHECK(Vector3<double>{1, 0,  0} == body.get_atom(2).coordinates());

        axis = {1, 1, 1};
        body.rotate(matrix::rotation_matrix(axis, std::numbers::pi/4));
        CHECK(Vector3<double>{-0.5058793634, 0.3106172175, -0.8047378541} == body.get_atom(0).coordinates());
        CHECK(Vector3<double>{-0.3106172175, 0.8047378541,  0.5058793634} == body.get_atom(1).coordinates());
        CHECK(Vector3<double>{ 0.8047378541, 0.5058793634, -0.3106172175} == body.get_atom(2).coordinates());
    }

    SECTION("simple 2") {
        std::vector<AtomFF> a = {
            AtomFF({1, 0, 0}, form_factor::form_factor_t::C),
            AtomFF({0, 1, 0}, form_factor::form_factor_t::C),
            AtomFF({0, 0, 1}, form_factor::form_factor_t::C)};
        Body body(a);

        body.rotate(matrix::rotation_matrix(0., std::numbers::pi/2, 0.));
        CHECK(Vector3<double>{0, 0, -1} == body.get_atom(0).coordinates());
        CHECK(Vector3<double>{0, 1,  0} == body.get_atom(1).coordinates());
        CHECK(Vector3<double>{1, 0,  0} == body.get_atom(2).coordinates());

        body.rotate(matrix::rotation_matrix(0.5612026, 0.3158423, 0.5612026));
        CHECK(Vector3<double>{-0.5058793634, 0.3106172175, -0.8047378541}.equals(body.get_atom(0).coordinates(), 1e-3));
        CHECK(Vector3<double>{-0.3106172175, 0.8047378541,  0.5058793634}.equals(body.get_atom(1).coordinates(), 1e-3));
        CHECK(Vector3<double>{ 0.8047378541, 0.5058793634, -0.3106172175}.equals(body.get_atom(2).coordinates(), 1e-3));
    }

    SECTION("complex") {
        std::vector<AtomFF> a = {
            AtomFF({0, 2, 1}, form_factor::form_factor_t::C),
            AtomFF({5, 1, 3}, form_factor::form_factor_t::C),
            AtomFF({6, 1, 4}, form_factor::form_factor_t::C),
            AtomFF({3, 7, 2}, form_factor::form_factor_t::C)};
        Body body(a);

        Vector3<double> axis = {0.5, 2, 1};
        body.rotate(matrix::rotation_matrix(axis, 1.8));
        REQUIRE(Vector3<double>{0.5843819499, 1.6706126346, 1.3665837559}  == body.get_atom(0).coordinates());
        REQUIRE(Vector3<double>{1.8656722055, 4.7666664324, -2.9661689675} == body.get_atom(1).coordinates());
        REQUIRE(Vector3<double>{2.6638285975, 5.6804357476, -3.692785794}  == body.get_atom(2).coordinates());
        REQUIRE(Vector3<double>{0.0886646879, 7.4409765368, 2.5737145825}  == body.get_atom(3).coordinates());
    }
}

TEST_CASE("Body::register_probe") {
    Body body(std::vector<AtomFF>{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)});
    state::StateManager sm(1);
    auto probe = std::make_shared<signaller::BoundSignaller>(0, &sm);
    body.register_probe(probe);
    REQUIRE(probe == body.get_signaller());
}

TEST_CASE_METHOD(multiple_fixture, "Body: equality") {
    SECTION("operator=") {
        b1 = b3;
        REQUIRE(b1.get_atoms().size() == 2);
        REQUIRE(b1.get_atom(0) == a5);
        REQUIRE(b1.get_atom(1) == a6);

        b1.get_atom(0) = a1;
        REQUIRE(b3.get_atom(0) == a5);

        // assignment with temporary bodies
        b1 = Body();
        {
            Body b5(std::vector<AtomFF>{a1, a2});
            b1 = b5;
        }
        REQUIRE(b1.get_atoms().size() == 2);
        REQUIRE(b1.get_atom(0) == a1);
        REQUIRE(b1.get_atom(1) == a2);
    }

    SECTION("equals_content") {
        b1 = b3;
        REQUIRE(b1.equals_content(b3));

        b2 = Body();
        {
            Body b5(std::vector<AtomFF>{a5, a6});
            b2 = b5;
        }
        REQUIRE(b2.equals_content(b3));

        b2.symmetry().add(symmetry::type::c3);
        REQUIRE(!b2.equals_content(b3));
        b4 = b3;
        REQUIRE(b3.equals_content(b4));
    }
}

TEST_CASE("Body::operator==") {
    std::vector<AtomFF> a1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> a2 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
    Body b1(a1);
    Body b2(a2);

    CHECK(!(b1 == b2)); // even though they have the same contents, body equality is defined exclusively by a uid
    Body b2c = b2;
    CHECK(b2 == b2c);
}

TEST_CASE_METHOD(fixture, "Body::state") {
    auto protein = Molecule({body});
    protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManager);
    auto* manager = static_cast<hist::IPartialHistogramManager*>(protein.get_histogram_manager())->get_state_manager();
    manager->reset_to_false();

    SECTION("Body::changed_external_state") {
        protein.get_body(0).get_signaller()->modified_external();
        CHECK(manager->get_externally_modified_bodies()[0] == true);
    }

    SECTION("Body::changed_internal_state") {
        protein.get_body(0).get_signaller()->modified_internal();
        CHECK(manager->get_internally_modified_bodies()[0] == true);
    }
}

TEST_CASE_METHOD(fixture, "Body::get_id") {
    int id = body.get_uid();
    Body b2(a);
    CHECK(id+1 == b2.get_uid());
}

TEST_CASE("Body::refresh_weights") {
    settings::general::verbose = false;
    SECTION("rebuilds the atom weights from the metadata") {
        settings::scattering::radiation = settings::scattering::Radiation::XRay;
        Body body("tests/files/2epe.pdb");
        auto xray = body.get_atoms();

        settings::scattering::radiation = settings::scattering::Radiation::Neutron;
        REQUIRE(body.refresh_weights());
        Body neutron("tests/files/2epe.pdb");
        REQUIRE(body.size_atom() == neutron.size_atom());
        for (int i = 0; i < body.size_atom(); ++i) {
            CHECK(body.get_atom(i).weight() == neutron.get_atom(i).weight());
            CHECK(body.get_atom(i).weight() != xray[i].weight());
        }
        settings::scattering::radiation = settings::scattering::Radiation::XRay;
    }

    SECTION("leaves weights it cannot derive untouched") {
        // a body built from bare atoms carries no metadata, so its weights are the user's and must not be guessed at
        Body body(std::vector<AtomFF>{AtomFF({0, 0, 0}, form_factor::form_factor_t::C, 42)});
        settings::scattering::radiation = settings::scattering::Radiation::Neutron;
        CHECK_FALSE(body.refresh_weights());
        CHECK(body.get_atom(0).weight() == 42);
        settings::scattering::radiation = settings::scattering::Radiation::XRay;
    }
}
