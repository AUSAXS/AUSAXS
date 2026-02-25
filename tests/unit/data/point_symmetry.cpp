// SPDX-License-Identifier: LGPL-3.0-or-later
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/PointSymmetry.h>

#include <cmath>

using namespace ausaxs;
using namespace ausaxs::symmetry;

TEST_CASE("PointSymmetry::constructor") {
    SECTION("default") {
        [[maybe_unused]] PointSymmetry s;
    }

    SECTION("translation and rotation") {
        PointSymmetry s({1, 2, 3}, {4, 5, 6});
        CHECK(s.translation == Vector3<double>{1, 2, 3});
        CHECK(s.rotation    == Vector3<double>{4, 5, 6});
    }
}

TEST_CASE("PointSymmetry::properties") {
    SECTION("always one repetition") {
        CHECK(PointSymmetry({0, 0, 0}, {0, 0, 1}).repetitions() == 1);
        CHECK(PointSymmetry({1, 2, 3}, {4, 5, 6}).repetitions() == 1);
    }

    SECTION("never closed") {
        CHECK_FALSE(PointSymmetry({0, 0, 0}, {0, 0, 0}).is_closed());
        CHECK_FALSE(PointSymmetry({1, 0, 0}, {0, 0, 1}).is_closed());
        CHECK_FALSE(PointSymmetry({0, 1, 0}, {0, 1, 0}).is_closed());
    }
}

TEST_CASE("PointSymmetry::clone") {
    PointSymmetry s({1, 2, 3}, {4, 5, 6});
    auto c = s.clone();
    auto* p = dynamic_cast<PointSymmetry*>(c.get());
    REQUIRE(p != nullptr);
    CHECK(p->translation == Vector3<double>{1, 2, 3});
    CHECK(p->rotation    == Vector3<double>{4, 5, 6});

    // clone is independent
    p->translation = {0, 0, 0};
    CHECK(s.translation == Vector3<double>{1, 2, 3});
}

TEST_CASE("PointSymmetry::add") {
    SECTION("translations add") {
        PointSymmetry s1({1, 0, 0}, {0, 0, 1});
        PointSymmetry s2({0, 2, 0}, {0, 1, 0});
        s1.add(&s2);
        CHECK(s1.translation == Vector3<double>{1, 2, 0});
        CHECK(s1.rotation    == Vector3<double>{0, 1, 1});
    }

    SECTION("zero translations") {
        PointSymmetry s1({0, 0, 0}, {1, 0, 0});
        PointSymmetry s2({0, 0, 0}, {-1, 0, 0});
        s1.add(&s2);
        CHECK(s1.translation == Vector3<double>{0, 0, 0});
        CHECK(s1.rotation    == Vector3<double>{0, 0, 0});
    }
}

TEST_CASE("PointSymmetry::span") {
    PointSymmetry s({1, 2, 3}, {4, 5, 6});
    auto st = s.span_translation();
    auto sr = s.span_rotation();

    REQUIRE(st.size() == 3);
    REQUIRE(sr.size() == 3);
    CHECK(st[0] == 1); CHECK(st[1] == 2); CHECK(st[2] == 3);
    CHECK(sr[0] == 4); CHECK(sr[1] == 5); CHECK(sr[2] == 6);

    // spans are live references – writes should update the underlying fields
    st[0] = 10;
    CHECK(s.translation.x() == 10);
    sr[2] = 99;
    CHECK(s.rotation.z() == 99);
}

TEST_CASE("PointSymmetry::get_transform") {
    const double c = std::cos(1.0), s = std::sin(1.0);
    const double tol = 1e-9;

    SECTION("pure translation of atom at cm") {
        // When v == cm, the rotation about cm leaves it fixed: copy = cm + translation
        PointSymmetry sym({3, 1, 0}, {0, 0, 0});
        auto f = sym.get_transform({2, 5, 0});
        auto v = f({2, 5, 0});
        CHECK(v == Vector3<double>(5, 6, 0));
    }

    SECTION("zero rotation vector gives identity") {
        // rotation={0,0,0} → R = identity → copy coincides with original (modulo translation)
        PointSymmetry sym({0, 0, 0}, {0, 0, 0});
        auto f = sym.get_transform({0, 0, 0});
        CHECK(f({1, 0, 0}) == Vector3<double>(1, 0, 0));
        CHECK(f({0, 1, 0}) == Vector3<double>(0, 1, 0));
        CHECK(f({0, 0, 1}) == Vector3<double>(0, 0, 1));
    }

    SECTION("rotation about z-axis (rotation={0,0,1})") {
        // rotation={0,0,gamma} → extrinsic Euler Rz(gamma)
        PointSymmetry sym({0, 0, 0}, {0, 0, 1});
        auto f = sym.get_transform({0, 0, 0});

        CHECK(f({1, 0, 0}) == Vector3<double>(c, s, 0));
        CHECK(f({0, 1, 0}) == Vector3<double>(-s, c, 0));
        CHECK(f({0, 0, 1}) == Vector3<double>(0, 0, 1));
    }

    SECTION("rotation about y-axis (rotation={0,1,0})") {
        // rotation={0,beta,0} → extrinsic Euler Ry(beta): (x,y,z)→(x*c+z*s, y, -x*s+z*c)
        PointSymmetry sym({0, 0, 0}, {0, 1, 0});
        auto f = sym.get_transform({0, 0, 0});

        CHECK(f({1, 0, 0}) == Vector3<double>( c, 0, -s));
        CHECK(f({0, 0, 1}) == Vector3<double>( s, 0,  c));
        CHECK(f({0, 1, 0}) == Vector3<double>(0, 1, 0));
    }

    SECTION("rotation about x-axis (rotation={1,0,0})") {
        // rotation={alpha,0,0} → extrinsic Euler Rx(alpha): (x,y,z)→(x, y*c-z*s, y*s+z*c)
        PointSymmetry sym({0, 0, 0}, {1, 0, 0});
        auto f = sym.get_transform({0, 0, 0});

        CHECK(f({0, 1, 0}) == Vector3<double>(0,  c, s));
        CHECK(f({0, 0, 1}) == Vector3<double>(0, -s, c));
        CHECK(f({1, 0, 0}) == Vector3<double>(1, 0, 0));
    }

    SECTION("rotation with non-zero cm: atom at cm maps to cm + translation") {
        // cm={2,0,0}, translation={0,0,0}: atom at cm stays at cm regardless of rotation
        PointSymmetry sym({0, 0, 0}, {0, 0, 1});
        auto f = sym.get_transform({2, 0, 0});
        CHECK(f({2, 0, 0}) == Vector3<double>(2, 0, 0));

        // atom at {3,0,0}: Rz(1)*(3-2) + 2 = Rz(1)*{1,0,0} + {2,0,0} = {c,s,0} + {2,0,0} = {2+c, s, 0}
        CHECK(f({3, 0, 0}) == Vector3<double>(2 + c, s, 0));
    }

    SECTION("combined translation and rotation") {
        // rotation={0,0,1}→Rz(1), translation={0,1,0}, cm={2,0,0}
        // f(cm) = cm + translation = {2,1,0}
        PointSymmetry sym({0, 1, 0}, {0, 0, 1});
        auto f = sym.get_transform({2, 0, 0});
        CHECK(f({2, 0, 0}) == Vector3<double>(2, 1, 0));

        // atom at {3,0,0}: Rz(1)*(3-2)+2+t = {c,s,0}+{2,1,0} = {2+c,1+s,0}
        // T = cm + t - R*cm = {2,0,0}+{0,1,0}-{2c,2s,0} = {2-2c, 1-2s, 0}
        // f({3,0,0}) = Rz(1)*{3,0,0} + T = {3c,3s,0} + {2-2c,1-2s,0} = {2+c, 1+s, 0}
        CHECK(f({3, 0, 0}) == Vector3<double>(2 + c, 1 + s, 0));
    }

    SECTION("different rotation magnitudes give different rotations") {
        // rotation = Euler angles: {0,0,1} → Rz(1 rad), {0,0,2} → Rz(2 rad)
        PointSymmetry s1({0, 0, 0}, {0, 0, 1});
        PointSymmetry s2({0, 0, 0}, {0, 0, 2});
        auto f1 = s1.get_transform({0, 0, 0});
        auto f2 = s2.get_transform({0, 0, 0});

        auto r1 = f1({1, 0, 0});
        auto r2 = f2({1, 0, 0});
        // Rz(1) and Rz(2) produce different results for the same input
        CHECK_FALSE((std::abs(r1.x() - r2.x()) < tol && std::abs(r1.y() - r2.y()) < tol));
    }

    SECTION("rep defaults to 1") {
        PointSymmetry sym({1, 0, 0}, {0, 0, 1});
        // rep=1 is the only valid value; non-default call must produce the same result
        auto f1 = sym.get_transform({0, 0, 0});
        auto f2 = sym.get_transform({0, 0, 0}, 1);
        auto r1 = f1({1, 2, 3});
        auto r2 = f2({1, 2, 3});
        CHECK(r1 == r2);
    }
}
