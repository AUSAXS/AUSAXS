// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/sequencer/detail/SymmetryFit.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/DihedralSymmetry.h>
#include <data/symmetry/PlanarDihedralSymmetry.h>
#include <data/symmetry/TetrahedralSymmetry.h>
#include <data/symmetry/OctahedralSymmetry.h>
#include <data/symmetry/IcosahedralSymmetry.h>
#include <data/symmetry/IPolyhedralSymmetry.h>
#include <math/MatrixUtils.h>

#include <numbers>
#include <random>

using namespace ausaxs;
using namespace ausaxs::symmetry;
using ausaxs::rigidbody::sequencer::detail::fit_symmetry;

namespace {
    std::vector<Vector3<double>> make_body(unsigned seed = 1, int n = 12) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(-5, 5);
        std::vector<Vector3<double>> v;
        for (int i = 0; i < n; ++i) {v.push_back({dist(gen), dist(gen), dist(gen)});}
        return v;
    }

    Vector3<double> centre_of(const std::vector<Vector3<double>>& body) {
        Vector3<double> cm{0, 0, 0};
        for (const auto& p : body) {cm += p;}
        return cm/static_cast<double>(body.size());
    }

    // expand a reference body through a source symmetry into the full set of copies (copies[0] = reference)
    std::vector<std::vector<Vector3<double>>> expand(const ISymmetry& source, const std::vector<Vector3<double>>& body) {
        auto cm = centre_of(body);
        std::vector<std::vector<Vector3<double>>> copies{body};
        for (unsigned int k = 1; k <= source.repetitions(); ++k) {
            auto t = source.get_transform(cm, k);
            std::vector<Vector3<double>> copy;
            for (const auto& p : body) {copy.push_back(t(p));}
            copies.push_back(std::move(copy));
        }
        return copies;
    }

    void add_noise(std::vector<std::vector<Vector3<double>>>& copies, double sigma, unsigned seed = 7) {
        std::mt19937 gen(seed);
        std::normal_distribution<double> nd(0, sigma);
        for (auto& c : copies) {for (auto& p : c) {p += Vector3<double>{nd(gen), nd(gen), nd(gen)};}}
    }
}

TEST_CASE("fit_symmetry::cyclic") {
    auto body = make_body();
    auto cm = centre_of(body);

    int n = GENERATE(2, 3, 4, 5, 6, 8, 12);
    Vector3<double> axis{0.3, 1.0, 0.5};
    Vector3<double> offset{4, 0, 0};
    CyclicSymmetry source(offset, {0, 0, 0}, axis, 2*std::numbers::pi/n, n - 1);

    auto copies = expand(source, body);
    REQUIRE(copies.size() == static_cast<std::size_t>(n));

    CyclicSymmetry templ({0, 0, 0}, {0, 0, 0}, {0, 0, 1}, 0, n - 1);
    auto res = fit_symmetry(templ, cm, copies);

    // the fitted symmetry reproduces the assembly
    CHECK_THAT(res.rmsd, Catch::Matchers::WithinAbs(0, 1e-6));

    // the recovered generator is a rotation by 2pi/n
    auto* cs = dynamic_cast<CyclicSymmetry*>(res.symmetry.get());
    REQUIRE(cs != nullptr);
    CHECK_THAT(cs->_repeat_relation.angle, Catch::Matchers::WithinAbs(2*std::numbers::pi/n, 1e-6));
}

TEST_CASE("fit_symmetry::point") {
    auto body = make_body();
    auto cm = centre_of(body);

    PointSymmetry source({6, 1, 2}, {0.4, 0.2, 1.1});
    auto copies = expand(source, body);

    PointSymmetry templ;
    auto res = fit_symmetry(templ, cm, copies);
    CHECK_THAT(res.rmsd, Catch::Matchers::WithinAbs(0, 1e-6));
}

TEST_CASE("fit_symmetry::dihedral") {
    auto body = make_body();
    auto cm = centre_of(body);

    auto check = [&](ISymmetry& source, const ISymmetry& templ) {
        auto copies = expand(source, body);
        auto res = fit_symmetry(templ, cm, copies);
        CHECK_THAT(res.rmsd, Catch::Matchers::WithinAbs(0, 1e-6));
    };

    SECTION("d2") {DihedralSymmetry<2> s, t; s.rotation = {0.5, -0.3, 0.9}; s.translation = {3, -2, 1}; check(s, t);}
    SECTION("d3") {DihedralSymmetry<3> s, t; s.rotation = {0.5, -0.3, 0.9}; s.translation = {3, -2, 1}; check(s, t);}
    SECTION("d4") {DihedralSymmetry<4> s, t; s.rotation = {0.5, -0.3, 0.9}; s.translation = {3, -2, 1}; check(s, t);}
    SECTION("d6") {DihedralSymmetry<6> s, t; s.rotation = {0.2,  0.8, -0.4}; s.translation = {1, 4, -2}; check(s, t);}
}

TEST_CASE("fit_symmetry::planar_dihedral") {
    auto body = make_body();
    auto cm = centre_of(body);

    auto check = [&](ISymmetry& source, const ISymmetry& templ) {
        auto copies = expand(source, body);
        auto res = fit_symmetry(templ, cm, copies);
        CHECK_THAT(res.rmsd, Catch::Matchers::WithinAbs(0, 1e-6));
    };

    // planar dihedral offsets live in the group plane (z-component 0)
    SECTION("dp3") {PlanarDihedralSymmetry<3> s, t; s.rotation = {0.5, -0.3, 0.9}; s.translation = {3, -2, 0}; check(s, t);}
    SECTION("dp4") {PlanarDihedralSymmetry<4> s, t; s.rotation = {0.5, -0.3, 0.9}; s.translation = {3, -2, 0}; check(s, t);}
    SECTION("dp6") {PlanarDihedralSymmetry<6> s, t; s.rotation = {0.1,  0.6, -0.7}; s.translation = {2,  5, 0}; check(s, t);}
}

TEST_CASE("fit_symmetry::polyhedral") {
    auto body = make_body();
    auto cm = centre_of(body);

    auto check = [&](IPolyhedralSymmetry& source, const ISymmetry& templ) {
        source.rotation = {0.5, -0.3, 0.9};
        source.translation = {3, -2, 1};
        auto copies = expand(source, body);
        auto res = fit_symmetry(templ, cm, copies);
        CHECK_THAT(res.rmsd, Catch::Matchers::WithinAbs(0, 1e-6));
    };

    SECTION("tetrahedral") {TetrahedralSymmetry s, t; check(s, t);}
    SECTION("octahedral")  {OctahedralSymmetry  s, t; check(s, t);}
    SECTION("icosahedral") {IcosahedralSymmetry s, t; check(s, t);}
}

TEST_CASE("fit_symmetry::degrades gracefully with noise") {
    auto body = make_body();
    auto cm = centre_of(body);
    double sigma = 0.05;

    SECTION("cyclic") {
        CyclicSymmetry source({4, 0, 0}, {0, 0, 0}, {0.3, 1, 0.5}, 2*std::numbers::pi/4, 3);
        auto copies = expand(source, body);
        add_noise(copies, sigma);
        CyclicSymmetry templ({0, 0, 0}, {0, 0, 0}, {0, 0, 1}, 0, 3);
        auto res = fit_symmetry(templ, cm, copies);
        // residual should stay on the order of the injected noise, not blow up
        CHECK(res.rmsd < 10*sigma);
    }

    SECTION("octahedral") {
        OctahedralSymmetry source, templ;
        source.rotation = {0.5, -0.3, 0.9};
        source.translation = {3, -2, 1};
        auto copies = expand(source, body);
        add_noise(copies, sigma);
        auto res = fit_symmetry(templ, cm, copies);
        CHECK(res.rmsd < 10*sigma);
    }
}
