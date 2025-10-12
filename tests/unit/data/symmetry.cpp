#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/Symmetry.h>

#include <numbers>

using namespace ausaxs;
using namespace ausaxs::symmetry;

TEST_CASE("Symmetry::Symmetry") {
    SECTION("default") {
        Symmetry s;
        CHECK(s.repetitions == 0);
    }

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
        Symmetry::_Relation r;
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
        Symmetry::_Repeat r;
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
    SECTION("translation only") {
        // Translate initial body by (1,2,3)
        Symmetry s({{1, 2, 3}, {0, 0, 0}});
        auto f = s.get_transform<double>({0, 0, 0});
        CHECK(f({1, 2, 3}) == Vector3<double>(2, 4, 6));
        CHECK(f({0, 0, 0}) == Vector3<double>(1, 2, 3));
        CHECK(f({-1, -1, -1}) == Vector3<double>(0, 1, 2));

        // Translate subsequent bodies by (0,1,0)
        s = Symmetry({{1, 0, 0}, {0, 0, 0}}, {{0, 1, 0}, {0, 0, 0}}, 4);
        f = s.get_transform<double>({0, 0, 0});
        CHECK(f({1, 2, 3}) == Vector3<double>(2, 3, 3));

        f = s.get_transform<double>({0, 0, 0}, 2);
        CHECK(f({1, 2, 3}) == Vector3<double>(2, 4, 3));

        f = s.get_transform<double>({0, 0, 0}, 3);
        CHECK(f({1, 2, 3}) == Vector3<double>(2, 5, 3));
    }

    SECTION("rotation about X by +90 deg") {
        // A rotation about X by +90° (π/2): 
        // (x, y, z) -> (x, z, -y)
        Symmetry s({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {std::numbers::pi/2, 0, 0}}, 1);
        auto f = s.get_transform<double>({0, 0, 0});

        CHECK(f({1, 0, 0}) == Vector3<double>(1,  0, 0));   // x-axis unchanged
        CHECK(f({0, 1, 0}) == Vector3<double>(0,  0, 1));   // y -> z
        CHECK(f({0, 0, 1}) == Vector3<double>(0, -1, 0));   // z -> -y
    }

    SECTION("rotation about Y by +90 deg") {
        // Rotate about Y by +90° (π/2):
        // (x, y, z) -> (z, y, -x)
        Symmetry s({{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, std::numbers::pi/2, 0}}, 1);
        auto f = s.get_transform<double>({0, 0, 0});

        CHECK(f({1, 0, 0}) == Vector3<double>(0, 0, -1));  // x -> -z
        CHECK(f({0, 1, 0}) == Vector3<double>(0, 1, 0));   // y unchanged
        CHECK(f({0, 0, 1}) == Vector3<double>(1, 0, 0));   // z -> x
    }

    SECTION("translation + rotation about X by +90 deg") {
        // Combine a translation (1,2,3) and internal rotation about X by +90°
        // Offset by (1,2,3)
        // Then rotate about X by +90°: (x, y, z) -> (x, z, -y)
        Symmetry s({{1, 2, 3}, {0, 0, 0}}, {{0, 0, 0}, {-std::numbers::pi/2, 0, 0}}, 1);

        {   // first copy
            auto f = s.get_transform<double>({0, 0, 0});

            // Translation: (1,0,0) -> (2,2,3)
            // Rotation:    (2,2,3) -> (2,3,-2)
            CHECK(f({1, 0, 0}) == Vector3<double>(2, 3, -2)); 

            // Translation: (0,1,0) -> (1,3,3)
            // Rotation:    (1,3,3) -> (1,3,-3)
            CHECK(f({0, 1, 0}) == Vector3<double>(1, 3, -3));

            // Translation: (0,0,1) -> (1,2,4)
            // Rotation:    (1,2,4) -> (1,4,-2)
            CHECK(f({0, 0, 1}) == Vector3<double>(1, 4, -2));
        }

        {   // second copy
            auto f = s.get_transform<double>({0, 0, 0}, 2);

            // Translation: (1,0,0) -> (2, 2, 3)
            // Rotation:    (2, 2, 3) -> (2, 3, -2) -> (2, -2, -3)
            CHECK(f({1, 0, 0}) == Vector3<double>(2, -2, -3));

            // Translation: (0,1,0) -> (1, 3, 3)
            // Rotation:    (1, 3, 3) -> (1, 3, -3) -> (1, -3, -3)
            CHECK(f({0, 1, 0}) == Vector3<double>(1, -3, -3));

            // Translation: (0,0,1) -> (1, 2, 4)
            // Rotation:    (1, 2, 4) -> (1, 4, -2) -> (1, -2, -4)
            CHECK(f({0, 0, 1}) == Vector3<double>(1, -2, -4));
        }
    }

    SECTION("rotation about Z by 180 deg + translation") {
        // Offset by (1, 1, 0)
        // Then offset by (0, 0, 1)
        // Then rotate about Z by 180°: (x, y, z) -> (-x, -y, z)
        Symmetry s({{1, 1, 0}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, -std::numbers::pi}}, 1);
        auto f = s.get_transform<double>({0, 0, 0});

        {   // first copy
            // Translation: (1,0,0) -> (2,1,1)
            // Rotation:    (2,1,1) -> (-2,-1,1)
            CHECK(f({1,0,0}) == Vector3<double>(-2,-1,1));

            // Translation: (0,1,0) -> (1,2,1)
            // Rotation:    (1,2,1) -> (-1,-2,1)
            CHECK(f({0,1,0}) == Vector3<double>(-1,-2,1));

            // Translation: (0,0,1) -> (1,1,2)
            // Rotation:    (1,1,1) -> (-1,-1,1)
            CHECK(f({0,0,1}) == Vector3<double>(-1,-1,2));
        }

        {   // second copy
            auto f = s.get_transform<double>({0, 0, 0}, 2);

            // Translation: (1,0,0) -> (2,1,2)
            // Rotation:    (2,1,2) -> (-2,-1,2) -> (2,1,2)
            CHECK(f({1,0,0}) == Vector3<double>(2,1,2));

            // Translation: (0,1,0) -> (1,2,2)
            // Rotation:    (1,2,2) -> (-1,-2,2) -> (1,2,2)
            CHECK(f({0,1,0}) == Vector3<double>(1,2,2));

            // Translation: (0,0,1) -> (1,1,3)
            // Rotation:    (1,1,2) -> (-1,-1,2) -> (1,1,2)
            CHECK(f({0,0,1}) == Vector3<double>(1,1,3));
        }
    }

    SECTION("identity") {
        Symmetry s({{0, 0, 0}, {0, 0, 0}});
        auto f = s.get_transform<double>({0, 0, 0});
        
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 2, 3));
        CHECK(f({0, 0, 0}) == Vector3<double>(0, 0, 0));
        CHECK(f({-1, -1, -1}) == Vector3<double>(-1, -1, -1));
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
