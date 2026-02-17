#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/Symmetry.h>

#include <numbers>

using namespace ausaxs;
using namespace ausaxs::symmetry;

TEST_CASE("Symmetry::Symmetry") {
    SECTION("_Relation") {
        Symmetry s({{1, 2, 3}, {0, 0, 0}});
        CHECK(s.initial_relation.translation == Vector3<double>{1, 2, 3});
        CHECK(s.initial_relation.orientation == Vector3<double>{0, 0, 0});
        CHECK(s.repetitions == 1);
    }

    SECTION("_Relation, _Repeat") {
        Symmetry s({{1, 0, 0}, {0, 0, 0}}, {{0, 1, 0}, {0, 0, 0}}, 5);
        CHECK(s.initial_relation.translation == Vector3<double>{1, 0, 0});
        CHECK(s.initial_relation.orientation == Vector3<double>{0, 0, 0});
        CHECK(s.repeat_relation.translate == Vector3<double>{0, 1, 0});
        CHECK(s.repeat_relation.rotate == Vector3<double>{0, 0, 0});
        CHECK(s.repetitions == 5);
    }

    SECTION("_Relation, repetitions") {
        Symmetry s({{0, 0, 1}, {std::numbers::pi/4, 0, 0}}, 3);
        CHECK(s.initial_relation.translation == Vector3<double>{0, 0, 1});
        CHECK_THAT(s.initial_relation.orientation.x(), Catch::Matchers::WithinAbs(std::numbers::pi/4, 1e-6));
        CHECK(s.repetitions == 3);
    }
}

TEST_CASE("Symmetry::_Relation") {
    SECTION("default") {
        [[maybe_unused]] Symmetry::_Relation r;
    }

    SECTION("rvalue construction") {
        Symmetry::_Relation r({1, 2, 3}, {4, 5, 6});
        CHECK(r.translation == Vector3<double>{1, 2, 3});
        CHECK(r.orientation == Vector3<double>{4, 5, 6});
    }

    SECTION("lvalue construction") {
        Vector3<double> t{1, 2, 3};
        Vector3<double> o{4, 5, 6};
        Symmetry::_Relation r(t, o);
        CHECK(r.translation == Vector3<double>{1, 2, 3});
        CHECK(r.orientation == Vector3<double>{4, 5, 6});
    }

    SECTION("equality") {
        Symmetry::_Relation r1({1, 2, 3}, {4, 5, 6});
        Symmetry::_Relation r2({1, 2, 3}, {4, 5, 6});
        Symmetry::_Relation r3({1, 2, 3}, {4, 5, 7});
        
        CHECK(r1 == r2);
        CHECK_FALSE(r1 == r3);
    }
}

TEST_CASE("Symmetry::_Repeat") {
    SECTION("default") {
        [[maybe_unused]] Symmetry::_Repeat r;
    }

    SECTION("rvalue construction") {
        Symmetry::_Repeat r({1, 2, 3}, {4, 5, 6});
        CHECK(r.translate == Vector3<double>{1, 2, 3});
        CHECK(r.rotate == Vector3<double>{4, 5, 6});
    }

    SECTION("lvalue construction") {
        Vector3<double> t{1, 2, 3};
        Vector3<double> r_vec{4, 5, 6};
        Symmetry::_Repeat r(t, r_vec);
        CHECK(r.translate == Vector3<double>{1, 2, 3});
        CHECK(r.rotate == Vector3<double>{4, 5, 6});
    }

    SECTION("equality") {
        Symmetry::_Repeat r1({1, 2, 3}, {4, 5, 6});
        Symmetry::_Repeat r2({1, 2, 3}, {4, 5, 6});
        Symmetry::_Repeat r3({1, 2, 4}, {4, 5, 6});
        
        CHECK(r1 == r2);
        CHECK_FALSE(r1 == r3);
    }
}

TEST_CASE("Symmetry::is_closed") {
    SECTION("translations") {
        SECTION("simple translations") {
            CHECK(Symmetry({{ 1,  0,  0}, {0, 0, 0}}).is_closed());
            CHECK(Symmetry({{ 0,  1,  0}, {0, 0, 0}}).is_closed());
            CHECK(Symmetry({{ 0,  0,  1}, {0, 0, 0}}).is_closed());
            CHECK(Symmetry({{ 1,  2,  3}, {0, 0, 0}}).is_closed());
            CHECK(Symmetry({{-1, -2, -3}, {0, 0, 0}}).is_closed());
        }

        SECTION("repeating translations") {
            CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{ 1,  0,  0}, {0, 0, 0}}, 5).is_closed());
            CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{ 0,  1,  0}, {0, 0, 0}}, 5).is_closed());
            CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{ 0,  0,  1}, {0, 0, 0}}, 5).is_closed());
            CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{ 1,  2,  3}, {0, 0, 0}}, 5).is_closed());
            CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{-1, -2, -3}, {0, 0, 0}}, 5).is_closed());
        }
    }

    SECTION("rotations") {
        SECTION("insufficient repeats") {
            SECTION("x") {
                int repeats = GENERATE(1, 2, 4, 5);
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/2, 0, 0}}, repeats).is_closed());
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/4, 0, 0}}, repeats).is_closed());
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/8, 0, 0}}, repeats).is_closed());
            }

            SECTION("y") {
                int repeats = GENERATE(1, 2, 4, 5);
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/2, 0}}, repeats).is_closed());
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/4, 0}}, repeats).is_closed());
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/8, 0}}, repeats).is_closed());
            }

            SECTION("z") {
                int repeats = GENERATE(1, 2, 4, 5);
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, std::numbers::pi/2}}, repeats).is_closed());
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, std::numbers::pi/4}}, repeats).is_closed());
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, std::numbers::pi/8}}, repeats).is_closed());
            }

            SECTION("sufficient repeats but translated") {
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {std::numbers::pi/2, 0, 0}}, 3).is_closed());
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 1, 0}, {0, std::numbers::pi/2, 0}}, 3).is_closed());
                CHECK_FALSE(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, std::numbers::pi/2}}, 3).is_closed());
            }
        }

        SECTION("sufficient repeats") {
            SECTION("x") {
                CHECK(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/2, 0, 0}}, 3).is_closed());
                CHECK(Symmetry({{1, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/3, 0, 0}}, 5).is_closed());
                CHECK(Symmetry({{0, 2, 0}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/4, 0, 0}}, 7).is_closed());
                CHECK(Symmetry({{0, 0, 3}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/5, 0, 0}}, 9).is_closed());
                CHECK(Symmetry({{1, 0, 3}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/6, 0, 0}}, 11).is_closed());
            }

            SECTION("y") {
                CHECK(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/2, 0}}, 3).is_closed());
                CHECK(Symmetry({{1, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/3, 0}}, 5).is_closed());
                CHECK(Symmetry({{0, 2, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/4, 0}}, 7).is_closed());
                CHECK(Symmetry({{0, 0, 3}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/5, 0}}, 9).is_closed());
                CHECK(Symmetry({{1, 0, 3}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/6, 0}}, 11).is_closed());
            }

            SECTION("z") {
                CHECK(Symmetry({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, std::numbers::pi/2}}, 3).is_closed());
                CHECK(Symmetry({{1, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, std::numbers::pi/3}}, 5).is_closed());
                CHECK(Symmetry({{0, 2, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, std::numbers::pi/4}}, 7).is_closed());
                CHECK(Symmetry({{0, 0, 3}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, std::numbers::pi/5}}, 9).is_closed());
                CHECK(Symmetry({{1, 0, 3}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, std::numbers::pi/6}}, 11).is_closed());
            }
        }
    }
}

TEST_CASE("Symmetry::get_transform") {
    SECTION("translation only, no rotation (degenerate)") {
        // With R_r = I, all copies overlap with the original (the rotation is the identity).
        // The transform should still work: T = R_r^k * t_i - t_i + sum = 0
        Symmetry s({{1, 2, 3}, {0, 0, 0}});
        auto f = s.get_transform<double>({0, 0, 0});
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 2, 3));
        CHECK(f({0, 0, 0}) == Vector3<double>(0, 0, 0));
        CHECK(f({-1, -1, -1}) == Vector3<double>(-1, -1, -1));

        // With a per-step translation t_r, copies are shifted by k*t_r from original
        s = Symmetry({{1, 0, 0}, {0, 0, 0}}, {{0, 1, 0}, {0, 0, 0}}, 4);
        f = s.get_transform<double>({0, 0, 0});
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 3, 3));

        f = s.get_transform<double>({0, 0, 0}, 2);
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 4, 3));

        f = s.get_transform<double>({0, 0, 0}, 3);
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 5, 3));
    }

    SECTION("rotation about X by +90 deg") {
        // With t_i = 0, copies rotate about the body's center (which is at the center of rotation)
        // Rotation about X by +90° (π/2): (x, y, z) -> (x, -z, y)
        Symmetry s({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/2, 0, 0}}, 1);
        auto f = s.get_transform<double>({0, 0, 0});

        CHECK(f({1, 0, 0}) == Vector3<double>(1,  0, 0));   // x-axis unchanged
        CHECK(f({0, 1, 0}) == Vector3<double>(0,  0, 1));   // y -> z
        CHECK(f({0, 0, 1}) == Vector3<double>(0, -1, 0));   // z -> -y
    }

    SECTION("rotation about Y by +90 deg") {
        // Rotate about Y by +90° (π/2): (x, y, z) -> (-z, y, x)
        Symmetry s({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/2, 0}}, 1);
        auto f = s.get_transform<double>({0, 0, 0});

        CHECK(f({1, 0, 0}) == Vector3<double>(0, 0, -1));  // x -> -z
        CHECK(f({0, 1, 0}) == Vector3<double>(0, 1, 0));   // y unchanged
        CHECK(f({0, 0, 1}) == Vector3<double>(1, 0, 0));   // z -> x
    }

    SECTION("offset + rotation about X by -90 deg") {
        // Body offset by t_i = (1,2,3), then copies generated by rotating -90° about X
        // R_r(-π/2, X): (x,y,z) -> (x, z, -y)
        // T = R_r*t_i - t_i = (1,3,-2) - (1,2,3) = (0, 1, -5)
        Symmetry s({{1, 2, 3}, {0, 0, 0}}, {{0, 0, 0}, {-std::numbers::pi/2, 0, 0}}, 1);

        {   // first copy (k=1)
            auto f = s.get_transform<double>({0, 0, 0});

            CHECK(f({1, 0, 0}) == Vector3<double>(1, 1, -5));
            CHECK(f({0, 1, 0}) == Vector3<double>(0, 1, -6));
            CHECK(f({0, 0, 1}) == Vector3<double>(0, 2, -5));
        }

        {   // second copy (k=2): R_r^2 = rotation(-π, X): (x,y,z) -> (x,-y,-z)
            // T_2 = R_r*T_1 + T_1 = R_r*(0,1,-5) + (0,1,-5) = (0,-5,-1) + (0,1,-5) = (0,-4,-6)
            auto f = s.get_transform<double>({0, 0, 0}, 2);

            CHECK(f({1, 0, 0}) == Vector3<double>(1, -4, -6));
            CHECK(f({0, 1, 0}) == Vector3<double>(0, -5, -6));
            CHECK(f({0, 0, 1}) == Vector3<double>(0, -4, -7));
        }
    }

    SECTION("rotation about Z by 180 deg + per-step translation") {
        // Body offset by t_i = (1,1,0), per-step translation t_r = (0,0,1)
        // R_r(-π, Z): (x,y,z) -> (-x,-y,z)
        // T_1 = R_r*t_i - t_i + t_r = (-1,-1,0) - (1,1,0) + (0,0,1) = (-2,-2,1)
        Symmetry s({{1, 1, 0}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, -std::numbers::pi}}, 1);

        {   // first copy (k=1)
            auto f = s.get_transform<double>({0, 0, 0});

            CHECK(f({1,0,0}) == Vector3<double>(-3,-2,1));
            CHECK(f({0,1,0}) == Vector3<double>(-2,-3,1));
            CHECK(f({0,0,1}) == Vector3<double>(-2,-2,2));
        }

        {   // second copy (k=2): R_r^2 = I for 180° rotation
            // T_2 = R_r*T_1 + T_base = (2,2,1) + (-2,-2,1) = (0,0,2)
            auto f = s.get_transform<double>({0, 0, 0}, 2);

            // f(v) = v + (0,0,2): copies shift along z by 2*t_r
            CHECK(f({1,0,0}) == Vector3<double>(1,0,2));
            CHECK(f({0,1,0}) == Vector3<double>(0,1,2));
            CHECK(f({0,0,1}) == Vector3<double>(0,0,3));
        }
    }

    SECTION("identity") {
        Symmetry s({{0, 0, 0}, {0, 0, 0}});
        auto f = s.get_transform<double>({0, 0, 0});
        
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 2, 3));
        CHECK(f({0, 0, 0}) == Vector3<double>(0, 0, 0));
        CHECK(f({-1, -1, -1}) == Vector3<double>(-1, -1, -1));
    }

    SECTION("p4 symmetry with non-zero CM") {
        // Body at cm=(10,20,0), offset by t_i=(5,0,0) from center of symmetry
        // Center of symmetry in real space: cm - t_i = (5, 20, 0)
        // All 4 members should be equidistant (distance 5) from center
        Vector3<double> cm{10, 20, 0};
        Symmetry s({{5, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, std::numbers::pi/2}}, 3);
        Vector3<double> center = cm - s.initial_relation.translation;  // (5, 20, 0)

        {   // Copy 1 (90° rotation)
            auto f = s.get_transform<double>(cm, 1);
            auto copy_cm = f(cm);
            CHECK(copy_cm == Vector3<double>(5, 25, 0));
            CHECK_THAT((copy_cm - center).magnitude(), Catch::Matchers::WithinAbs(5, 1e-6));
        }

        {   // Copy 2 (180° rotation)
            auto f = s.get_transform<double>(cm, 2);
            auto copy_cm = f(cm);
            CHECK(copy_cm == Vector3<double>(0, 20, 0));
            CHECK_THAT((copy_cm - center).magnitude(), Catch::Matchers::WithinAbs(5, 1e-6));
        }

        {   // Copy 3 (270° rotation)
            auto f = s.get_transform<double>(cm, 3);
            auto copy_cm = f(cm);
            CHECK(copy_cm == Vector3<double>(5, 15, 0));
            CHECK_THAT((copy_cm - center).magnitude(), Catch::Matchers::WithinAbs(5, 1e-6));
        }

        // Original at cm=(10,20,0) is also at distance 5 from center (5,20,0)
        CHECK_THAT((cm - center).magnitude(), Catch::Matchers::WithinAbs(5, 1e-6));
    }
}

TEST_CASE("Symmetry::equality") {
    Symmetry s1({{1, 2, 3}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}, 5);
    Symmetry s2({{1, 2, 3}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}, 5);
    Symmetry s3({{1, 2, 4}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}, 5);
    Symmetry s4({{1, 2, 3}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}, 6);
    
    CHECK(s1 == s2);
    CHECK_FALSE(s1 == s3);
    CHECK_FALSE(s1 == s4);
}
