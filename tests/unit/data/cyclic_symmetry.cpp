#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/CyclicSymmetry.h>

#include <numbers>

using namespace ausaxs;
using namespace ausaxs::symmetry;

TEST_CASE("Symmetry::Symmetry") {
    SECTION("_Relation, _Repeat") {
        // t_r = {0,0,1} lies along the rotation axis {0,0,1} (screw symmetry)
        CyclicSymmetry s({{1, 0, 0}}, {{0, 0, 1}, {0, 0, 1}, 0.5}, 5);
        CHECK(s._initial_relation.translation == Vector3<double>{1, 0, 0});
        CHECK(s._repeat_relation.translation == Vector3<double>{0, 0, 1});
        CHECK(s._repeat_relation.axis        == Vector3<double>{0, 0, 1});
        CHECK_THAT(s._repeat_relation.angle, Catch::Matchers::WithinAbs(0.5, 1e-9));
        CHECK(s._repetitions == 5);
    }
}

TEST_CASE("Symmetry::_Relation") {
    SECTION("default") {
        [[maybe_unused]] CyclicSymmetry::_Relation r;
    }

    SECTION("construction") {
        CyclicSymmetry::_Relation r({1, 2, 3});
        CHECK(r.translation == Vector3<double>{1, 2, 3});
    }

    SECTION("equality") {
        CyclicSymmetry::_Relation r1({1, 2, 3});
        CyclicSymmetry::_Relation r2({1, 2, 3});
        CyclicSymmetry::_Relation r3({1, 2, 4});

        CHECK(r1 == r2);
        CHECK_FALSE(r1 == r3);
    }
}

TEST_CASE("Symmetry::_Repeat") {
    SECTION("default") {
        [[maybe_unused]] CyclicSymmetry::_Repeat r;
    }

    SECTION("full construction") {
        CyclicSymmetry::_Repeat r({1, 0, 0}, {0, 0, 1}, 1.5);
        CHECK(r.translation == Vector3<double>{1, 0, 0});
        CHECK(r.axis        == Vector3<double>{0, 0, 1});
        CHECK_THAT(r.angle, Catch::Matchers::WithinAbs(1.5, 1e-9));
    }

    SECTION("axis+angle constructor (no translation)") {
        CyclicSymmetry::_Repeat r({0, 1, 0}, std::numbers::pi);
        CHECK(r.translation == Vector3<double>{0, 0, 0});
        CHECK(r.axis        == Vector3<double>{0, 1, 0});
        CHECK_THAT(r.angle, Catch::Matchers::WithinAbs(std::numbers::pi, 1e-9));
    }

    SECTION("equality") {
        CyclicSymmetry::_Repeat r1({1, 0, 0}, {0, 0, 1}, 1.5);
        CyclicSymmetry::_Repeat r2({1, 0, 0}, {0, 0, 1}, 1.5);
        CyclicSymmetry::_Repeat r3({1, 0, 0}, {0, 0, 1}, 1.6);

        CHECK(r1 == r2);
        CHECK_FALSE(r1 == r3);
    }
}

TEST_CASE("Symmetry::is_closed") {
    SECTION("translations") {
        SECTION("simple translations (no rotation)") {
            CHECK(CyclicSymmetry({ 1,  0,  0}, {0, 0, 0}, {0, 0, 1}, 0).is_closed());
            CHECK(CyclicSymmetry({ 0,  1,  0}, {0, 0, 0}, {0, 0, 1}, 0).is_closed());
            CHECK(CyclicSymmetry({ 0,  0,  1}, {0, 0, 0}, {0, 0, 1}, 0).is_closed());
            CHECK(CyclicSymmetry({ 1,  2,  3}, {0, 0, 0}, {0, 0, 1}, 0).is_closed());
            CHECK(CyclicSymmetry({-1, -2, -3}, {0, 0, 0}, {0, 0, 1}, 0).is_closed());
        }

        SECTION("repeating translations block closure") {
            CHECK_FALSE(CyclicSymmetry({{0,0,0}}, {{1,0,0}, {0,1,0}, 0}, 5).is_closed());
            CHECK_FALSE(CyclicSymmetry({{0,0,0}}, {{0,1,0}, {0,1,0}, 0}, 5).is_closed());
            CHECK_FALSE(CyclicSymmetry({{0,0,0}}, {{0,0,1}, {0,0,1}, 0}, 5).is_closed());
        }
    }

    SECTION("rotations") {
        SECTION("insufficient repeats") {
            int repeats = GENERATE(1, 2, 4, 5);
            CHECK_FALSE(CyclicSymmetry({{0,0,0}}, {{1,0,0}, std::numbers::pi/2}, repeats).is_closed());
            CHECK_FALSE(CyclicSymmetry({{0,0,0}}, {{0,1,0}, std::numbers::pi/2}, repeats).is_closed());
            CHECK_FALSE(CyclicSymmetry({{0,0,0}}, {{0,0,1}, std::numbers::pi/2}, repeats).is_closed());
        }

        SECTION("sufficient repeats but with translation") {
            CHECK_FALSE(CyclicSymmetry({{0,0,0}}, {{1, 0, 0}, {1, 0, 0}, std::numbers::pi/2}, 3).is_closed());
            CHECK_FALSE(CyclicSymmetry({{0,0,0}}, {{0, 1, 0}, {0, 1, 0}, std::numbers::pi/2}, 3).is_closed());
            CHECK_FALSE(CyclicSymmetry({{0,0,0}}, {{0, 0, 1}, {0, 0, 1}, std::numbers::pi/2}, 3).is_closed());
        }

        SECTION("sufficient repeats") {
            SECTION("x-axis") {
                CHECK(CyclicSymmetry({{0,0,0}}, {{1,0,0}, std::numbers::pi/2}, 3).is_closed());
                CHECK(CyclicSymmetry({{1,0,0}}, {{1,0,0}, std::numbers::pi/3}, 5).is_closed());
                CHECK(CyclicSymmetry({{0,2,0}}, {{1,0,0}, std::numbers::pi/4}, 7).is_closed());
            }

            SECTION("y-axis") {
                CHECK(CyclicSymmetry({{0,0,0}}, {{0,1,0}, std::numbers::pi/2}, 3).is_closed());
                CHECK(CyclicSymmetry({{1,0,0}}, {{0,1,0}, std::numbers::pi/3}, 5).is_closed());
            }

            SECTION("z-axis") {
                CHECK(CyclicSymmetry({{0,0,0}}, {{0,0,1}, std::numbers::pi/2}, 3).is_closed());
                CHECK(CyclicSymmetry({{1,0,0}}, {{0,0,1}, std::numbers::pi/3}, 5).is_closed());
                CHECK(CyclicSymmetry({{0,2,0}}, {{0,0,1}, std::numbers::pi/4}, 7).is_closed());
                CHECK(CyclicSymmetry({{0,0,3}}, {{0,0,1}, std::numbers::pi/5}, 9).is_closed());
                CHECK(CyclicSymmetry({{1,0,3}}, {{0,0,1}, std::numbers::pi/6}, 11).is_closed());
            }

            SECTION("predefined symmetries are closed") {
                CHECK(CyclicSymmetry({{0,0,0}}, {{0,0,1}, std::numbers::pi},   1).is_closed()); // p2
                CHECK(CyclicSymmetry({{0,0,0}}, {{0,0,1}, std::numbers::pi*2./3}, 2).is_closed()); // p3
                CHECK(CyclicSymmetry({{0,0,0}}, {{0,0,1}, std::numbers::pi/2}, 3).is_closed()); // p4
            }
        }
    }
}

TEST_CASE("Symmetry::get_transform") {
    SECTION("translation only, no rotation (degenerate)") {
        // Per-step translation shifts copies along y
        CyclicSymmetry s({1, 2, 3}, {0, 0, 0}, {0, 0, 0}, 0, 4);
        s = CyclicSymmetry({{1, 0, 0}}, {{0, 1, 0}, {0, 0, 1}, 0}, 4);
        auto f = s.get_transform({0, 0, 0});
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 3, 3));

        f = s.get_transform({0, 0, 0}, 2);
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 4, 3));

        f = s.get_transform({0, 0, 0}, 3);
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 5, 3));
    }

    SECTION("rotation about X by +90 deg") {
        // R_x(+pi/2): (x,y,z) -> (x,-z,y)
        CyclicSymmetry s({{0,0,0}}, {{1,0,0}, std::numbers::pi/2}, 1);
        auto f = s.get_transform({0, 0, 0});

        CHECK(f({1, 0, 0}) == Vector3<double>(1,  0, 0));
        CHECK(f({0, 1, 0}) == Vector3<double>(0,  0, 1));
        CHECK(f({0, 0, 1}) == Vector3<double>(0, -1, 0));
    }

    SECTION("rotation about Y by +90 deg") {
        // R_y(+pi/2): (x,y,z) -> (z,y,-x)
        CyclicSymmetry s({{0,0,0}}, {{0,1,0}, std::numbers::pi/2}, 1);
        auto f = s.get_transform({0, 0, 0});

        CHECK(f({1, 0, 0}) == Vector3<double>(0, 0, -1));
        CHECK(f({0, 1, 0}) == Vector3<double>(0, 1,  0));
        CHECK(f({0, 0, 1}) == Vector3<double>(1, 0,  0));
    }

    SECTION("offset + rotation about X by -90 deg") {
        // Body offset t_i=(1,2,3), R_r=R_x(-pi/2): (x,y,z)->(x,z,-y)
        // T_1 = R_r*t_i - t_i = (1,3,-2)-(1,2,3) = (0,1,-5)
        CyclicSymmetry s({{1, 2, 3}}, {{1, 0, 0}, -std::numbers::pi/2}, 1);

        {   // first copy (k=1)
            auto f = s.get_transform({0, 0, 0});
            CHECK(f({1, 0, 0}) == Vector3<double>(1, 1, -5));
            CHECK(f({0, 1, 0}) == Vector3<double>(0, 1, -6));
            CHECK(f({0, 0, 1}) == Vector3<double>(0, 2, -5));
        }

        {   // second copy (k=2): R_r^2 = R_x(-pi) -> (x,y,z)->(x,-y,-z)
            // T_2 = R_r*T_1 + T_1 = (0,-5,-1) + (0,1,-5) = (0,-4,-6)
            auto f = s.get_transform({0, 0, 0}, 2);
            CHECK(f({1, 0, 0}) == Vector3<double>(1, -4, -6));
            CHECK(f({0, 1, 0}) == Vector3<double>(0, -5, -6));
            CHECK(f({0, 0, 1}) == Vector3<double>(0, -4, -7));
        }
    }

    SECTION("rotation about Z by 180 deg + per-step translation") {
        // t_i=(1,1,0), t_r=(0,0,1), R_z(-pi): (x,y,z)->(-x,-y,z)
        // T_1 = R_r*t_i - t_i + t_r = (-1,-1,0)-(1,1,0)+(0,0,1) = (-2,-2,1)
        CyclicSymmetry s({{1, 1, 0}}, {{0, 0, 1}, {0, 0, 1}, -std::numbers::pi}, 1);

        {   // first copy (k=1)
            auto f = s.get_transform({0, 0, 0});
            CHECK(f({1,0,0}) == Vector3<double>(-3,-2,1));
            CHECK(f({0,1,0}) == Vector3<double>(-2,-3,1));
            CHECK(f({0,0,1}) == Vector3<double>(-2,-2,2));
        }

        {   // second copy (k=2): R_r^2=I, T_2=R_r*T_1+T_base=(2,2,1)+(-2,-2,1)=(0,0,2)
            auto f = s.get_transform({0, 0, 0}, 2);
            CHECK(f({1,0,0}) == Vector3<double>(1,0,2));
            CHECK(f({0,1,0}) == Vector3<double>(0,1,2));
            CHECK(f({0,0,1}) == Vector3<double>(0,0,3));
        }
    }

    SECTION("p4 symmetry with non-zero CM") {
        // Body at cm=(10,20,0), t_i=(5,0,0). Center of symmetry = (5,20,0).
        // All 4 members equidistant (r=5) from center.
        Vector3<double> cm{10, 20, 0};
        CyclicSymmetry s({{5, 0, 0}}, {{0, 0, 1}, std::numbers::pi/2}, 3);
        Vector3<double> center = cm - s._initial_relation.translation;  // (5,20,0)

        auto check_copy = [&](int rep, Vector3<double> expected) {
            auto f = s.get_transform(cm, rep);
            auto copy_cm = f(cm);
            CHECK(copy_cm == expected);
            CHECK_THAT((copy_cm - center).magnitude(), Catch::Matchers::WithinAbs(5, 1e-6));
        };

        check_copy(1, {5, 25, 0});
        check_copy(2, {0, 20, 0});
        check_copy(3, {5, 15, 0});
        CHECK_THAT((cm - center).magnitude(), Catch::Matchers::WithinAbs(5, 1e-6));
    }

    SECTION("rotation axis direction: non-z axes") {
        Vector3<double> cm{0, 0, 0};

        // p2 around X: (x,y,z) -> (x,-y,-z)
        CyclicSymmetry s_x({{0,0,0}}, {{1,0,0}, std::numbers::pi}, 1);
        auto p = s_x.get_transform(cm)({3, 4, 5});
        CHECK_THAT(p.x(), Catch::Matchers::WithinAbs( 3, 1e-9));
        CHECK_THAT(p.y(), Catch::Matchers::WithinAbs(-4, 1e-9));
        CHECK_THAT(p.z(), Catch::Matchers::WithinAbs(-5, 1e-9));

        // p2 around Y: (x,y,z) -> (-x,y,-z)
        CyclicSymmetry s_y({{0,0,0}}, {{0,1,0}, std::numbers::pi}, 1);
        auto q = s_y.get_transform(cm)({3, 4, 5});
        CHECK_THAT(q.x(), Catch::Matchers::WithinAbs(-3, 1e-9));
        CHECK_THAT(q.y(), Catch::Matchers::WithinAbs( 4, 1e-9));
        CHECK_THAT(q.z(), Catch::Matchers::WithinAbs(-5, 1e-9));

        // Axis magnitude should not matter: {7,0,0} same as {1,0,0}
        CyclicSymmetry s_x2({{0,0,0}}, {{7,0,0}, std::numbers::pi}, 1);
        auto r = s_x2.get_transform(cm)({3, 4, 5});
        CHECK_THAT(r.x(), Catch::Matchers::WithinAbs(p.x(), 1e-9));
        CHECK_THAT(r.y(), Catch::Matchers::WithinAbs(p.y(), 1e-9));
        CHECK_THAT(r.z(), Catch::Matchers::WithinAbs(p.z(), 1e-9));
    }
}