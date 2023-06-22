#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <utility/Constants.h>
#include <utility/Console.h>
#include <hydrate/Grid.h>
#include <hydrate/GridMember.h>
#include <data/Protein.h>
#include <data/state/StateManager.h>
#include <data/BodySplitter.h>
#include <settings/All.h>
#include <data/Water.h>
#include <data/Body.h>
#include <hist/HistogramManager.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using std::cout, std::endl, std::vector, std::shared_ptr;

struct fixture {
    vector<Atom> a = {Atom(1, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "0"), Atom(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, -1), 1, 0, "C", "0"),
                      Atom(3, "C", "", "LYS", "", 1, "", Vector3<double>( 1, -1, -1), 1, 0, "C", "0"), Atom(4, "C", "", "LYS", "", 1, "", Vector3<double>( 1, 1, -1), 1, 0, "C", "0"),
                      Atom(5, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1,  1), 1, 0, "C", "0"), Atom(6, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1,  1), 1, 0, "C", "0"),
                      Atom(7, "C", "", "LYS", "", 1, "", Vector3<double>( 1, -1,  1), 1, 0, "C", "0"), Atom(8, "C", "", "LYS", "", 1, "", Vector3<double>( 1, 1,  1), 1, 0, "C", "0")
    };
    Body body = Body(a, {});
};

struct multiple_fixture {
    Atom a1 = Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1);
    Atom a2 = Atom(Vector3<double>(-1,  1, -1), 1, "C", "C", 1);
    Atom a3 = Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1);
    Atom a4 = Atom(Vector3<double>(-1,  1,  1), 1, "C", "C", 1);
    Atom a5 = Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1);
    Atom a6 = Atom(Vector3<double>( 1,  1, -1), 1, "C", "C", 1);
    Atom a7 = Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1);
    Atom a8 = Atom(Vector3<double>( 1,  1,  1), 1, "He", "He", 1);
    Water w1 = Water(Vector3<double>(0, 1, 2), 1, "H", "H", 1);
    Water w2 = Water(Vector3<double>(3, 4, 5), 1, "H", "H", 1);

    Body b1 = Body(std::vector<Atom>{a1, a2});
    Body b2 = Body(std::vector<Atom>{a3, a4});
    Body b3 = Body(std::vector<Atom>{a5, a6});
    Body b4 = Body(std::vector<Atom>{a7, a8});
    std::vector<Body> ap = {b1, b2, b3, b4};
    Protein protein = Protein(ap);
};

TEST_CASE_METHOD(multiple_fixture, "Body::Body") {
    SECTION("ExistingFile&") {
        io::ExistingFile ef("test/files/2epe.pdb");
        Body b(ef);
        REQUIRE(b.get_atoms().size() == 260);
    }

    SECTION("vector<Atom>&") {
        Body b(std::vector<Atom>{a1, a2});
        REQUIRE(b.get_atoms().size() == 2);
        CHECK(b.get_atom(0) == a1);
        CHECK(b.get_atom(1) == a2);
    }

    SECTION("vector<Atom>&, vector<Water>&") {
        Body b(std::vector<Atom>{a1, a2}, std::vector<Water>{w1, w2});
        REQUIRE(b.get_atoms().size() == 2);
        CHECK(b.get_atom(0) == a1);
        CHECK(b.get_atom(1) == a2);
        REQUIRE(b.get_waters().size() == 2);
        CHECK(b.get_waters()[0] == w1);
        CHECK(b.get_waters()[1] == w2);
    }

    SECTION("Body&") {
        Body b(b1);
        REQUIRE(b.get_atoms().size() == 2);
        CHECK(b.get_atom(0) == a1);
        CHECK(b.get_atom(1) == a2);

        // check that they are backed by separate files
        b.get_atom(0) = a3;
        CHECK(b1.get_atom(0) == a1);
    }

    SECTION("Body&&") {
        Body b5 = std::move(b1);
        REQUIRE(b5.get_atoms().size() == 2);
        CHECK(b5.get_atom(0) == a1);
        CHECK(b5.get_atom(1) == a2);
    }
}

TEST_CASE("Body::save") {
    vector<Atom> a = {Atom(1, "C"  , "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, -1), 1, 0, "C", "0"),
                      Atom(3, "O"  , "", "LYS", "", 1, "", Vector3<double>(1, -1, -1), 1, 0, "O", "0"),   Atom(4, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, -1), 1, 0, "C", "0"),
                      Atom(5, "N"  , "", "LYS", "", 1, "", Vector3<double>(-1, -1, 1), 1, 0, "N", "0"),   Atom(6, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, 1), 1, 0, "C", "0"),
                      Atom(7, "OXT", "", "LYS", "", 1, "", Vector3<double>(1, -1, 1), 1, 0, "O", "0"),    Atom(8, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, 1), 1, 0, "C", "0")};
    Body body(a, {});

    body.save("temp/body_io.pdb");
    Body body2("temp/body_io.pdb");

    CHECK(body.get_atoms().size() == body2.get_atoms().size());
    for (unsigned int i = 0; i < body.get_atoms().size(); i++) {
        CHECK(body.get_atom(i).as_pdb() == body2.get_atom(i).as_pdb());
    }

    remove("temp/body_io.pdb");
}

TEST_CASE_METHOD(fixture, "Body::get_atoms") {
    CHECK(body.get_atoms() == a);
}

TEST_CASE_METHOD(fixture, "Body::get_atom") {
    REQUIRE(body.get_atoms().size() == a.size());
    for (unsigned int i = 0; i < body.get_atoms().size(); i++) {
        CHECK(body.get_atom(i) == a[i]);
    }
}

TEST_CASE_METHOD(multiple_fixture, "Body::get_waters") {
    auto waters = std::vector<Water>{w1, w2};
    Body body(std::vector<Atom>{a1, a2}, waters);
    REQUIRE(body.get_waters() == waters);
}

TEST_CASE_METHOD(fixture, "Body::get_cm") {
    Vector3<double> cm = body.get_cm();
    REQUIRE(cm == Vector3<double>({0, 0, 0}));
}

TEST_CASE("Body::get_volume_acids") {
    CHECK(false);
}

TEST_CASE("Body::get_volume_calpha") {
    CHECK(false);
}

TEST_CASE_METHOD(fixture, "Body::get_volume") {
    REQUIRE(body.get_volume_acids() == constants::volume::amino_acids.get("LYS"));
}


TEST_CASE_METHOD(fixture, "Body::molar_mass") {
    CHECK_THAT(body.molar_mass(), Catch::Matchers::WithinRel(8*constants::mass::atomic.get("C")*constants::Avogadro, 1e-6));
    console::print_warning("Check definition of molar mass.");
}

TEST_CASE_METHOD(fixture, "Body::absolute_mass") {
    CHECK_THAT(body.absolute_mass(), Catch::Matchers::WithinRel(8*constants::mass::atomic.get("C"), 1e-6));
}

TEST_CASE_METHOD(fixture, "Body::total_atomic_charge") {
    CHECK(body.total_atomic_charge() == 8*6);
}

TEST_CASE_METHOD(fixture, "Body::total_effective_charge") {
    double c0 = body.get_atom(0).get_effective_charge();    
    double charge1 = body.total_effective_charge();
    body.update_effective_charge(1.5);
    double charge2 = body.total_effective_charge();
    CHECK(charge2 == charge1+12);
    CHECK(body.get_atom(0).get_effective_charge() == c0+1.5);
}

TEST_CASE_METHOD(fixture, "Body::center") {
    SECTION("trivial center") {
        body.translate(Vector3<double>{-1, -1, -1});
        body.center();
        CHECK(body.get_atom(0).coords == Vector3<double>{-1, -1, -1});
        CHECK(body.get_atom(1).coords == Vector3<double>{-1,  1, -1});
        CHECK(body.get_atom(2).coords == Vector3<double>{ 1, -1, -1});
        CHECK(body.get_atom(3).coords == Vector3<double>{ 1,  1, -1});
    }

    SECTION("non-trivial center") {
        auto a = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "O", "O", 1),
             Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "O", "O", 1)};
        auto body = Body(a, {});

        body.center();
        double shift = 0.142402;

        // make checks in a roundabout fashion for better failure messages
        auto res = Vector3<double>{-1, -1-shift, 0};
        if (body.get_atom(0).coords.equals(res, 1e-3)) {CHECK(true);}
        else {CHECK(body.get_atom(0).coords == res);}

        res = Vector3<double>{-1, 1-shift, 0};
        if (body.get_atom(1).coords.equals(res, 1e-3)) {CHECK(true);}
        else {CHECK(body.get_atom(1).coords == res);}

        res = Vector3<double>{1, -1-shift, 0};
        if (body.get_atom(2).coords.equals(res, 1e-3)) {CHECK(true);}
        else {CHECK(body.get_atom(2).coords == res);}

        res = Vector3<double>{1, 1-shift, 0};
        if (body.get_atom(3).coords.equals(res, 1e-3)) {CHECK(true);}
        else {CHECK(body.get_atom(3).coords == res);}
    }
}

TEST_CASE_METHOD(fixture, "Body::translate") {
    SECTION("basic translation") {
        body.translate(Vector3<double>{1, 1, 1});
        CHECK(body.get_atom(0).coords == Vector3<double>{0, 0, 0});
        CHECK(body.get_atom(1).coords == Vector3<double>{0, 2, 0});
        CHECK(body.get_atom(2).coords == Vector3<double>{2, 0, 0});
        CHECK(body.get_atom(3).coords == Vector3<double>{2, 2, 0});
    }

    SECTION("informs manager") {
        auto protein = Protein({body});
        auto manager = protein.get_histogram_manager()->get_state_manager();
        manager->reset();
        protein.get_body(0).translate(Vector3<double>(10, 0, 0));
        CHECK(protein.get_body(0).get_atom(0).coords == Vector3<double>(9, -1, -1));
        CHECK(protein.get_body(0).get_atom(1).coords == Vector3<double>(9, 1, -1));
        CHECK(manager->get_externally_modified_bodies()[0] == true);
    }
}

TEST_CASE("Body::rotate") {
    SECTION("Matrix<double>&") {}

    SECTION("Vector3<double>&, double") {
        SECTION("simple") {
            vector<Atom> a = {Atom(Vector3<double>(1, 0, 0), 1, "C", "C", 1), 
                            Atom(Vector3<double>(0, 1, 0), 1, "C", "C", 1), 
                            Atom(Vector3<double>(0, 0, 1), 1, "C", "C", 1)};
            Body body(a, {});

            Vector3<double> axis = {0, 1, 0};
            body.rotate(axis, M_PI_2);
            CHECK(Vector3<double>({0, 0, -1}) == body.get_atom(0).coords); 
            CHECK(Vector3<double>({0, 1, 0}) == body.get_atom(1).coords); 
            CHECK(Vector3<double>({1, 0, 0}) == body.get_atom(2).coords); 

            axis = {1, 1, 1};
            body.rotate(axis, M_PI/4);
            CHECK(Vector3<double>({-0.5058793634, 0.3106172175, -0.8047378541}) == body.get_atom(0).coords); 
            CHECK(Vector3<double>({-0.3106172175, 0.8047378541, 0.5058793634}) == body.get_atom(1).coords); 
            CHECK(Vector3<double>({0.8047378541, 0.5058793634, -0.3106172175}) == body.get_atom(2).coords); 
        }

        SECTION("complex") {
            vector<Atom> a = {Atom(Vector3<double>(0, 2, 1), 1, "C", "C", 1), 
                            Atom(Vector3<double>(5, 1, 3), 1, "C", "C", 1), 
                            Atom(Vector3<double>(6, 1, 4), 1, "C", "C", 1),
                            Atom(Vector3<double>(3, 7, 2), 1, "C", "C", 1)};
            Body body(a, {});

            Vector3<double> axis = {0.5, 2, 1};
            body.rotate(axis, 1.8);
            REQUIRE(Vector3<double>({0.5843819499, 1.6706126346, 1.3665837559}) == body.get_atom(0).coords); 
            REQUIRE(Vector3<double>({1.8656722055, 4.7666664324, -2.9661689675}) == body.get_atom(1).coords); 
            REQUIRE(Vector3<double>({2.6638285975, 5.6804357476, -3.692785794}) == body.get_atom(2).coords); 
            REQUIRE(Vector3<double>({0.0886646879, 7.4409765368, 2.5737145825}) == body.get_atom(3).coords); 
        }
    }

    SECTION("double, double, double, ") {
        SECTION("simple") {
            vector<Atom> a = {Atom(Vector3<double>(1, 0, 0), 1, "C", "C", 1), 
                            Atom(Vector3<double>(0, 1, 0), 1, "C", "C", 1), 
                            Atom(Vector3<double>(0, 0, 1), 1, "C", "C", 1)};
            Body body(a, {});

            body.rotate(0, M_PI_2, 0);
            CHECK(Vector3<double>({0, 0, -1}) == body.get_atom(0).coords); 
            CHECK(Vector3<double>({0, 1, 0}) == body.get_atom(1).coords); 
            CHECK(Vector3<double>({1, 0, 0}) == body.get_atom(2).coords); 

            body.rotate(0.5612026, 0.3158423, 0.5612026);
            CHECK(Vector3<double>({-0.5058793634, 0.3106172175, -0.8047378541}).equals(body.get_atom(0).coords, 1e-3));
            CHECK(Vector3<double>({-0.3106172175, 0.8047378541, 0.5058793634}).equals(body.get_atom(1).coords, 1e-3));
            CHECK(Vector3<double>({0.8047378541, 0.5058793634, -0.3106172175}).equals(body.get_atom(2).coords, 1e-3));
        }
    }
}

TEST_CASE("Body::register_probe") {
    CHECK(false);
}

TEST_CASE_METHOD(multiple_fixture, "Body::operator=") {
    b1 = b3;
    REQUIRE(b1.get_atoms().size() == 2);
    REQUIRE(b1.get_atom(0) == a5);
    REQUIRE(b1.get_atom(1) == a6);

    b1.get_atom(0) = a1;
    REQUIRE(b3.get_atom(0) == a5);

    // assignment with temporary bodies
    b1 = Body();
    {
        Body b5({a1, a2});
        b1 = b5;
    }
    REQUIRE(b1.get_atoms().size() == 2);
    REQUIRE(b1.get_atom(0) == a1);
    REQUIRE(b1.get_atom(1) == a2);

    CHECK(false);
}

TEST_CASE("Body::operator==") {
    std::vector<Atom> a1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    std::vector<Atom> a2 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    Body b1(a1);
    Body b2(a2);

    CHECK(!(b1 == b2)); // even though they have the same contents, body equality is defined exclusively by a uid
    Body b2c = b2;
    CHECK(b2 == b2c);
}

TEST_CASE("Body::get_file") {
    CHECK(false);
}

TEST_CASE_METHOD(fixture, "Body::state") {
    auto protein = Protein({body});
    auto manager = protein.get_histogram_manager()->get_state_manager();
    manager->reset();

    SECTION("Body::changed_external_state") {
        protein.get_body(0).changed_external_state();
        CHECK(manager->get_externally_modified_bodies()[0] == true);
    }

    SECTION("Body::changed_internal_state") {
        protein.get_body(0).changed_internal_state();
        CHECK(manager->get_internally_modified_bodies()[0] == true);
    }
}

TEST_CASE_METHOD(fixture, "Body::get_id") {
    unsigned int id = body.get_id();
    Body b2;
    CHECK(id+1 == b2.get_id());
}

TEST_CASE_METHOD(fixture, "Body::atom_size") {
    CHECK(body.atom_size() == 8);
    CHECK(Body().atom_size() == 0);
}

TEST_CASE("grid") {
    SECTION("single") {
        Body b({Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)});
        grid::Grid g(Axis3D(-2, 2, -2, 2, -2, 2, 1));

        g.add(&b);
        REQUIRE(g.a_members.size() == 2);
        REQUIRE(g.get_volume() != 0);

        g.remove(&b);
        REQUIRE(g.a_members.size() == 0);
        REQUIRE(g.get_volume() == 0);
    }

    SECTION("multiple") {
        settings::grid::ra = 1;
        settings::grid::rh = 1;
        vector<Atom> a1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
        vector<Atom> a2 = {Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
        vector<Atom> a3 = {Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
        vector<Atom> a4 = {Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
        Body b1(a1), b2(a2), b3(a3), b4(a4);
        vector<Body> bodies = {b1, b2, b3, b4};
        grid::Grid grid(Axis3D(-5, 5, -5, 5, -5, 5, 1));

        grid.add(&b1);
        grid.add(&b2);
        grid.add(&b3);
        grid.add(&b4);
        REQUIRE(grid.a_members.size() == 8);

        unsigned int vol = grid.get_volume();
        grid.remove(&b2);
        grid.add(&b2);
        REQUIRE(grid.get_volume() == vol);
        grid.remove(&b2);
        grid.force_expand_volume();
        REQUIRE(grid.a_members.size() == 6);

        vol = grid.get_volume();
        grid.remove(&b1);
        grid.add(&b1);
        REQUIRE(grid.get_volume() == vol);
        grid.remove(&b1);
        grid.force_expand_volume();
        REQUIRE(grid.a_members.size() == 4);

        vol = grid.get_volume();
        grid.remove(&b3);
        grid.add(&b3);
        REQUIRE(grid.get_volume() == vol);
        grid.remove(&b3);
        grid.force_expand_volume();
        REQUIRE(grid.a_members.size() == 2);

        auto remaining = grid.a_members;
        for (const auto& e : remaining) {
            REQUIRE((e == a4[0] || e == a4[1]));
        }

        // check volume
        REQUIRE(grid.get_volume() != 0);
        grid.remove(&b4);
        REQUIRE(grid.get_volume() == 0);
    }

    SECTION("real data") {
        settings::grid::ra = 1;
        settings::grid::rh = 1;
        Protein protein = rigidbody::BodySplitter::split("test/files/2epe.pdb", {9, 99});
        unsigned int N = protein.get_atoms().size();
        auto grid = protein.get_grid();
        REQUIRE(grid->a_members.size() == N);
        REQUIRE(grid->get_volume() != 0);

        unsigned int vol = grid->get_volume();
        grid->remove(&protein.get_body(0));
        grid->add(&protein.get_body(0));
        REQUIRE(grid->get_volume() == vol);
        grid->remove(&protein.get_body(0));
        grid->force_expand_volume();
        REQUIRE(grid->a_members.size() == N - protein.get_body(0).get_atoms().size());
        REQUIRE(grid->get_volume() != 0);

        vol = grid->get_volume();
        grid->remove(&protein.get_body(1));
        grid->add(&protein.get_body(1));
        REQUIRE(grid->get_volume() == vol);
        grid->remove(&protein.get_body(1));        
        grid->force_expand_volume();
        REQUIRE(grid->a_members.size() == N - protein.get_body(0).get_atoms().size() - protein.get_body(1).get_atoms().size());
        REQUIRE(grid->get_volume() != 0);

        vol = grid->get_volume();
        grid->remove(&protein.get_body(2));
        grid->add(&protein.get_body(2));
        REQUIRE(grid->get_volume() == vol);
        grid->remove(&protein.get_body(2));
        REQUIRE(grid->a_members.size() == 0);
        REQUIRE(grid->get_volume() == 0);
    }

    // SECTION("perfect reset") {
    //     vector<Atom> a1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    //     vector<Atom> a2 = {Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
    //     vector<Atom> a3 = {Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
    //     vector<Atom> a4 = {Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
    //     Body b1(a1), b2(a2), b3(a3), b4(a4);
    //     vector<Body> bodies = {b1, b2, b3, b4};
    //     grid::Grid grid(Axis3D(-50, 50, -50, 50, -50, 50, 1));

    //     REQUIRE(grid.get_volume() == 0);
    //     grid.add(&b1);
    //     grid.add(&b2);
    //     grid.add(&b3);
    //     grid.add(&b4);
    //     REQUIRE(grid.a_members.size() == 8);

    //     auto grid_copy = grid;
    //     grid.remove(&b1);
    //     grid.remove(&b2);
    //     b1.translate(Vector3<double>(0, 0, 10));
    //     b2.translate(Vector3<double>(0, 0, 10));
    //     grid.add(&b1);
    //     grid.add(&b2);
    //     REQUIRE(grid.a_members.size() == 8);
    //     REQUIRE(grid.get_volume() == 0);

    //     // do the reset
    //     grid = grid_copy;
    //     REQUIRE(grid.a_members.size() == 8);

    //     grid.add(&b1);
    //     grid.add(&b2);
    //     grid.add(&b3);
    //     grid.add(&b4);
    //     REQUIRE(grid.a_members.size() == 8);
    //     REQUIRE(grid.get_volume() != 0);
    // }
}