#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/CompositeDistanceHistogram.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <io/ExistingFile.h>
#include <utility/Utility.h>

hist::CompositeDistanceHistogram generate_random(unsigned int size) {
    std::vector<double> p_pp(size), p_hp(size), p_hh(size), p(size);
    for (unsigned int i = 0; i < size; ++i) {
        p_pp[i] = rand() % 100;
        p_hp[i] = rand() % 100;
        p_hh[i] = rand() % 100;
        p[i] = p_pp[i] + p_hp[i] + p_hh[i];
    }
    Axis axis(1, 10, 1);
    return hist::CompositeDistanceHistogram(std::move(p_pp), std::move(p_hp), std::move(p_hh), std::move(p), axis);
}

TEST_CASE("CompositeDistanceHistogram::reset_water_scaling_factor") {
    auto hist = generate_random(100);
    auto p = hist.p;
    hist.apply_water_scaling_factor(2);
    CHECK(hist.p != p);
    hist.reset_water_scaling_factor();
    CHECK(hist.p == p);
}

TEST_CASE("CompositeDistanceHistogram::apply_water_scaling_factor") {
    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    std::vector<Atom> b1 =   {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1),  Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    std::vector<Atom> b2 =   {Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1),  Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
    std::vector<Atom> b3 =   {Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1),  Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
    std::vector<Water> w =   {Water(Vector3<double>(1, -1,  1), 1, "C", "C", 1),  Water(Vector3<double>(1, 1,  1), 1, "C", "C", 1)};
    std::vector<Body> a = {Body(b1), Body(b2), Body(b3)};
    Protein protein(a, w);

    auto hist = protein.get_histogram();
    std::vector<double> p_pp = hist->get_pp_histogram();
    std::vector<double> p_hp = hist->get_hp_histogram();
    std::vector<double> p_hh = hist->get_hh_histogram();

    hist->apply_water_scaling_factor(2);
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + 2*p_hp[i] + 4*p_hh[i], Catch::Matchers::WithinRel(hist->p[i]));
    }

    hist->apply_water_scaling_factor(3);
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + 3*p_hp[i] + 9*p_hh[i], Catch::Matchers::WithinRel(hist->p[i]));
    }

    hist->reset_water_scaling_factor();
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + p_hp[i] + p_hh[i], Catch::Matchers::WithinRel(hist->p[i]));
    }
}