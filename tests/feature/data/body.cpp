#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <constants/Constants.h>
#include <utility/Console.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <data/Molecule.h>
#include <data/state/StateManager.h>
#include <rigidbody/BodySplitter.h>
#include <settings/All.h>
#include <data/Body.h>
#include <hist/histogram_manager/IPartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <io/Writer.h>

#include <vector>
#include <string>
#include <numbers>

using namespace ausaxs;
using namespace data;

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
    multiple_fixture() {
        settings::molecule::center = false;
    }

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

TEST_CASE_METHOD(multiple_fixture, "Body::Body") {
    settings::general::verbose = false;
    SECTION("ExistingFile&") {
        io::ExistingFile ef("tests/files/2epe.pdb");
        Body b(ef);
        REQUIRE(b.size_atom() == 1001);
        REQUIRE(b.size_water() == 48);
    }

    SECTION("vector<Atom>&") {
        Body b(std::vector<AtomFF>{a1, a2});
        REQUIRE(b.size_atom() == 2);
        CHECK(b.get_atom(0) == a1);
        CHECK(b.get_atom(1) == a2);
    }

    SECTION("vector<Atom>&, vector<Water>&") {
        Body b(std::vector<AtomFF>{a1, a2}, std::vector<Water>{w1, w2});
        REQUIRE(b.size_atom() == 2);
        CHECK(b.get_atom(0) == a1);
        CHECK(b.get_atom(1) == a2);
        REQUIRE(b.size_water() == 2);
        CHECK(b.get_waters()->get()[0] == w1);
        CHECK(b.get_waters()->get()[1] == w2);
    }

    SECTION("Body&") {
        Body b(b1);
        REQUIRE(b.size_atom() == 2);
        CHECK(b.get_atom(0) == a1);
        CHECK(b.get_atom(1) == a2);

        // check that they are backed by separate files
        b.get_atom(0) = a3;
        CHECK(b1.get_atom(0) == a1);
    }

    SECTION("Body&&") {
        Body b5 = std::move(b1);
        REQUIRE(b5.size_atom() == 2);
        CHECK(b5.get_atom(0) == a1);
        CHECK(b5.get_atom(1) == a2);
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
    io::Writer::write({body}, "temp/body_io.pdb");
    Body body2("temp/body_io.pdb");

    CHECK(body.size_atom() == body2.size_atom());
    for (unsigned int i = 0; i < body.size_atom(); i++) {
        CHECK(body.get_atom(i) == body2.get_atom(i));
    }

    std::remove("temp/body_io.pdb");
}

TEST_CASE_METHOD(fixture, "Body::translate") {
    SECTION("basic translation") {
        body.translate(Vector3<double>{1, 1, 1});
        CHECK(body.get_atom(0).coordinates() == Vector3<double>{0, 0, 0});
        CHECK(body.get_atom(1).coordinates() == Vector3<double>{0, 2, 0});
        CHECK(body.get_atom(2).coordinates() == Vector3<double>{2, 0, 0});
        CHECK(body.get_atom(3).coordinates() == Vector3<double>{2, 2, 0});
    }

    SECTION("informs manager") {
        auto protein = Molecule({body});
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManager);
        auto manager = static_cast<hist::IPartialHistogramManager*>(protein.get_histogram_manager())->get_state_manager();
        manager->reset_to_false();
        protein.get_body(0).translate(Vector3<double>(10, 0, 0));
        CHECK(protein.get_body(0).get_atom(0).coordinates() == Vector3<double>(9, -1, -1));
        CHECK(protein.get_body(0).get_atom(1).coordinates() == Vector3<double>(9,  1, -1));
        CHECK(manager->get_externally_modified_bodies()[0] == true);
    }
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
        CHECK(Vector3<double>({0, 0, -1}) == body.get_atom(0).coordinates()); 
        CHECK(Vector3<double>({0, 1,  0}) == body.get_atom(1).coordinates()); 
        CHECK(Vector3<double>({1, 0,  0}) == body.get_atom(2).coordinates()); 

        axis = {1, 1, 1};
        body.rotate(matrix::rotation_matrix(axis, std::numbers::pi/4));
        CHECK(Vector3<double>({-0.5058793634, 0.3106172175, -0.8047378541}) == body.get_atom(0).coordinates()); 
        CHECK(Vector3<double>({-0.3106172175, 0.8047378541,  0.5058793634}) == body.get_atom(1).coordinates()); 
        CHECK(Vector3<double>({ 0.8047378541, 0.5058793634, -0.3106172175}) == body.get_atom(2).coordinates()); 
    }

    SECTION("simple 2") {
        std::vector<AtomFF> a = {
            AtomFF({1, 0, 0}, form_factor::form_factor_t::C), 
            AtomFF({0, 1, 0}, form_factor::form_factor_t::C), 
            AtomFF({0, 0, 1}, form_factor::form_factor_t::C)};
        Body body(a);

        body.rotate(matrix::rotation_matrix(0., std::numbers::pi/2, 0.));
        CHECK(Vector3<double>({0, 0, -1}) == body.get_atom(0).coordinates()); 
        CHECK(Vector3<double>({0, 1,  0}) == body.get_atom(1).coordinates()); 
        CHECK(Vector3<double>({1, 0,  0}) == body.get_atom(2).coordinates()); 

        body.rotate(matrix::rotation_matrix(0.5612026, 0.3158423, 0.5612026));
        CHECK(Vector3<double>({-0.5058793634, 0.3106172175, -0.8047378541}).equals(body.get_atom(0).coordinates(), 1e-3));
        CHECK(Vector3<double>({-0.3106172175, 0.8047378541,  0.5058793634}).equals(body.get_atom(1).coordinates(), 1e-3));
        CHECK(Vector3<double>({ 0.8047378541, 0.5058793634, -0.3106172175}).equals(body.get_atom(2).coordinates(), 1e-3));
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
        REQUIRE(Vector3<double>({0.5843819499, 1.6706126346, 1.3665837559})  == body.get_atom(0).coordinates()); 
        REQUIRE(Vector3<double>({1.8656722055, 4.7666664324, -2.9661689675}) == body.get_atom(1).coordinates()); 
        REQUIRE(Vector3<double>({2.6638285975, 5.6804357476, -3.692785794})  == body.get_atom(2).coordinates()); 
        REQUIRE(Vector3<double>({0.0886646879, 7.4409765368, 2.5737145825})  == body.get_atom(3).coordinates()); 
    }
}

#include <data/state/BoundSignaller.h>
TEST_CASE("Body::register_probe") {
    Body body(std::vector<AtomFF>{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)});
    auto probe = std::make_shared<signaller::BoundSignaller>(1, nullptr);
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

        b2.symmetry().add(symmetry::type::p3);
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
    auto manager = static_cast<hist::IPartialHistogramManager*>(protein.get_histogram_manager())->get_state_manager();
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

TEST_CASE_METHOD(fixture, "Body::atom_size") {
    CHECK(body.size_atom() == 8);
    CHECK(Body().size_atom() == 0);
}