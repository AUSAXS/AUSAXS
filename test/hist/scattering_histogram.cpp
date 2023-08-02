#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/ScatteringHistogram.h>

ScatteringHistogram generate_random(unsigned int size) {
    std::vector<double> p_pp(size), p_hp(size), p_hh(size), p(size);
    for (unsigned int i = 0; i < size; ++i) {
        p_pp[i] = rand() % 100;
        p_hp[i] = rand() % 100;
        p_hh[i] = rand() % 100;
        p[i] = p_pp[i] + p_hp[i] + p_hh[i];
    }
    Axis axis(1, 10, 1);
    return ScatteringHistogram(p_pp, p_hp, p_hh, p, axis);
}

TEST_CASE("ScatteringHistogram::ScatteringHistogram") {
    SECTION("default") {
        hist::ScatteringHistogram hist;
        CHECK(hist.size() == 0);
        CHECK(hist.limits() == Limit{0, 0});
    }

    SECTION("ScatteringHistogram&&") {
        auto hist = generate_random(10);
        auto hist_old = hist;
        auto hist_new = std::move(hist);
        CHECK(hist_new == hist_old);
    }

    SECTION("ScatteringHistogram&") {
        auto hist = generate_random(10);
        auto hist_new = hist;
        CHECK(hist_new == hist);
    }

    SECTION("vector<double>&, vector<double>&, vector<double>&, vector<double>&, Axis&") {
        std::vector<double> p_pp{1, 2, 3, 4, 5}, p_hp{1, 2, 3, 4, 5}, p_hh{1, 2, 3, 4, 5}, p{1, 2, 3, 4, 5};
        Axis axis(1, 10, 1);
        hist::ScatteringHistogram hist(p_pp, p_hp, p_hh, p, axis);
        CHECK(hist.size() == 5);
        CHECK(hist.limits() == Limit{1, 10});
    }
}

TEST_CASE("ScatteringHistograM::calc_guinier_gyration_ratio_squared") {CHECK(false);}

TEST_CASE("ScatteringHistogram::calc_debye_scattering_intensity") {
    SECTION("default") {CHECK(false);}
    SECTION("vector<double>&") {CHECK(false);}
}

TEST_CASE("ScatteringHistogram::operator=") {
    SECTION("copy") {
        auto hist = generate_random(10);
        auto hist_new = hist;
        CHECK(hist_new == hist);

    }

    SECTION("move") {
        auto hist = generate_random(10);
        auto hist_old = hist;
        auto hist_new = std::move(hist);
        CHECK(hist_new == hist_old);
    }
}

TEST_CASE("ScatteringHistogram::operator+=") {
    auto hist1 = generate_random(100);
    auto hist2 = generate_random(100);
    auto hist1_old = hist1;
    auto hist2_old = hist2;
    hist1 += hist2;
    for (unsigned int i = 0; i < hist1.size(); i++) {
        CHECK(hist1.p_pp[i] == hist1_old.p_pp[i] + hist2_old.p_pp[i]);
        CHECK(hist1.p_hp[i] == hist1_old.p_hp[i] + hist2_old.p_hp[i]);
        CHECK(hist1.p_hh[i] == hist1_old.p_hh[i] + hist2_old.p_hh[i]);
        CHECK(hist1.p[i] == hist1_old.p[i] + hist2_old.p[i]);
    }
}

TEST_CASE("ScatteringHistogram::operator-=") {
    auto hist1 = generate_random(100);
    auto hist2 = generate_random(100);
    auto hist1_old = hist1;
    auto hist2_old = hist2;
    hist1 -= hist2;
    for (unsigned int i = 0; i < hist1.size(); i++) {
        CHECK(hist1.p_pp[i] == hist1_old.p_pp[i] - hist2_old.p_pp[i]);
        CHECK(hist1.p_hp[i] == hist1_old.p_hp[i] - hist2_old.p_hp[i]);
        CHECK(hist1.p_hh[i] == hist1_old.p_hh[i] - hist2_old.p_hh[i]);
        CHECK(hist1.p[i] == hist1_old.p[i] - hist2_old.p[i]);
    }    
}

TEST_CASE("ScatteringHistogram::operator*=") {
    auto hist1 = generate_random(100);
    auto hist1_old = hist1;
    hist1 *= 2;
    for (unsigned int i = 0; i < hist1.size(); i++) {
        CHECK(hist1.p_pp[i] == hist1_old.p_pp[i] * 2);
        CHECK(hist1.p_hp[i] == hist1_old.p_hp[i] * 2);
        CHECK(hist1.p_hh[i] == hist1_old.p_hh[i] * 2);
        CHECK(hist1.p[i] == hist1_old.p[i] * 2);
    }
}

TEST_CASE("ScatteringHistogram::operator==") {
    auto hist1 = generate_random(100);
    auto hist2 = hist1;
    CHECK(hist1 == hist2);

    hist2 = generate_random(100);
    CHECK(hist1 != hist2);
}

TEST_CASE("ScatteringHistogram::extend_axis") {CHECK(false);}

TEST_CASE("ScatteringHistogram::reset_water_scaling_factor") {
    auto hist = generate_random(100);
    auto p = hist.p;
    hist.apply_water_scaling_factor(2);
    CHECK(hist.p != p);
    hist.reset_water_scaling_factor();
    CHECK(hist.p == p);
}

TEST_CASE("ScatteringHistogram::apply_water_scaling_factor") {
    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    std::vector<Atom> b1 =   {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1),  Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    std::vector<Atom> b2 =   {Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1),  Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
    std::vector<Atom> b3 =   {Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1),  Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
    std::vector<Water> w =   {Water(Vector3<double>(1, -1,  1), 1, "C", "C", 1),  Water(Vector3<double>(1, 1,  1), 1, "C", "C", 1)};
    std::vector<Body> a = {Body(b1), Body(b2), Body(b3)};
    Protein protein(a, w);

    hist::ScatteringHistogram hist = protein.get_histogram();
    std::vector<double> p_pp = hist.p_pp.p;
    std::vector<double> p_hp = hist.p_hp.p;
    std::vector<double> p_hh = hist.p_hh.p;

    hist.apply_water_scaling_factor(2);
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + 2*p_hp[i] + 4*p_hh[i], Catch::Matchers::WithinRel(hist.p[i]));
    }

    hist.apply_water_scaling_factor(3);
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + 3*p_hp[i] + 9*p_hh[i], Catch::Matchers::WithinRel(hist.p[i]));
    }

    hist.reset_water_scaling_factor();
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + p_hp[i] + p_hh[i], Catch::Matchers::WithinRel(hist.p[i]));
    }
}