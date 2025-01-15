#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Symmetry.h>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("Symmetry::is_closed") {
    SECTION("translation only") {
        auto s = GENERATE(
            data::detail::Symmetry{{1, 0, 0}},
            data::detail::Symmetry{{0, 1, 0}},
            data::detail::Symmetry{{0, 0, 1}},
            data::detail::Symmetry{{1, 2, 3}},
            data::detail::Symmetry{{-1, -2, -3}}
        );
        CHECK_FALSE(s.is_closed()); // translations are never closed
    }

    SECTION("rotations") {
        SECTION("insufficient repeats") {
            SECTION("x") {
                int repeats = GENERATE(1, 2, 4, 5);
                CHECK_FALSE(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {std::numbers::pi/2, 0, 0}, repeats).is_closed());
                CHECK_FALSE(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {std::numbers::pi/4, 0, 0}, repeats).is_closed());
                CHECK_FALSE(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {std::numbers::pi/8, 0, 0}, repeats).is_closed());
            }

            SECTION("y") {
                int repeats = GENERATE(1, 2, 4, 5);
                CHECK_FALSE(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/2, 0}, repeats).is_closed());
                CHECK_FALSE(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/4, 0}, repeats).is_closed());
                CHECK_FALSE(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/8, 0}, repeats).is_closed());
            }

            SECTION("z") {
                int repeats = GENERATE(1, 2, 4, 5);
                CHECK_FALSE(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/2}, repeats).is_closed());
                CHECK_FALSE(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/4}, repeats).is_closed());
                CHECK_FALSE(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/8}, repeats).is_closed());
            }

            SECTION("sufficient repeats but translated") {
                CHECK_FALSE(data::detail::Symmetry({1, 0, 0}, {0, 0, 0}, {std::numbers::pi/2, 0, 0}, 3).is_closed());
                CHECK_FALSE(data::detail::Symmetry({1, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/2, 0}, 3).is_closed());
                CHECK_FALSE(data::detail::Symmetry({1, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/2}, 3).is_closed());
            }
        }

        SECTION("sufficient repeats") {
            SECTION("x") {
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {std::numbers::pi/2, 0, 0}, 3).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {std::numbers::pi/3, 0, 0}, 5).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {std::numbers::pi/4, 0, 0}, 7).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {std::numbers::pi/5, 0, 0}, 9).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {std::numbers::pi/6, 0, 0}, 11).is_closed());
            }

            SECTION("y") {
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/2, 0}, 3).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/3, 0}, 5).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/4, 0}, 7).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/5, 0}, 9).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/6, 0}, 11).is_closed());
            }

            SECTION("z") {
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/2}, 3).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/3}, 5).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/4}, 7).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/5}, 9).is_closed());
                CHECK(data::detail::Symmetry({0, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/6}, 11).is_closed());
            }
        }
    }
}

TEST_CASE("Symmetry::get_transform") {
    SECTION("no transformation") {
        // Default constructor: no translation, no rotation
        data::detail::Symmetry s;
        auto f = s.get_transform<double>();
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 2, 3));
    }

    SECTION("translation only") {
        // Translate by (1,2,3)
        data::detail::Symmetry s({1, 2, 3});
        auto f = s.get_transform<double>();
        CHECK(f({1, 2, 3}) == Vector3<double>(2, 4, 6));
        CHECK(f({0, 0, 0}) == Vector3<double>(1, 2, 3));
        CHECK(f({-1, -1, -1}) == Vector3<double>(0, 1, 2));
    }

    SECTION("rotation about X by +90 deg") {
        // A rotation about X by +90° (π/2): 
        // (x, y, z) -> (x, z, -y)
        data::detail::Symmetry s({0, 0, 0}, {0, 0, 0}, {std::numbers::pi/2, 0, 0}, 1);
        auto f = s.get_transform<double>();
        
        CHECK(f({1, 0, 0}) == Vector3<double>(1,  0, 0));   // x-axis unchanged
        CHECK(f({0, 1, 0}) == Vector3<double>(0,  0, 1));   // y -> z
        CHECK(f({0, 0, 1}) == Vector3<double>(0, -1, 0));   // z -> -y
    }

    SECTION("rotation about Y by +90 deg") {
        // Rotate about Y by +90° (π/2):
        // (x, y, z) -> (z, y, -x)
        data::detail::Symmetry s({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/2, 0}, 1);
        auto f = s.get_transform<double>();
        
        CHECK(f({1, 0, 0}) == Vector3<double>(0, 0, -1));  // x -> -z
        CHECK(f({0, 1, 0}) == Vector3<double>(0, 1, 0));   // y unchanged
        CHECK(f({0, 0, 1}) == Vector3<double>(1, 0, 0));   // z -> x
    }

    SECTION("translation + rotation about X by +90 deg") {
        // Combine a translation (1,2,3) and internal rotation about X by +90°
        // Rotation about X by +90°: (x, y, z) -> (x, z, -y)
        // Then translate by (1,2,3).
        data::detail::Symmetry s({1, 2, 3}, {0, 0, 0}, {std::numbers::pi/2, 0, 0}, 1);
        auto f = s.get_transform<double>();
        
        // Before translation: (1,0,0) -> (1,0,0)
        // After translation: (2,2,3)
        CHECK(f({1, 0, 0}) == Vector3<double>(2, 2, 3));

        // Before translation: (0,1,0) -> (0,0,1)
        // After translation: (1,2,4)
        CHECK(f({0, 1, 0}) == Vector3<double>(1, 2, 4));

        // Before translation: (0,0,1) -> (0,-1,0)
        // After translation: (1,1,3)
        CHECK(f({0, 0, 1}) == Vector3<double>(1, 1, 3));
    }

    SECTION("rotation about Z by 180 deg + translation") {
        // Rotate about Z by 180°: (x, y, z) -> (-x, -y, z)
        // Then translate by (1,1,0).
        data::detail::Symmetry s({1, 1, 0}, {0, 0, 0}, {0, 0, std::numbers::pi}, 1);
        auto f = s.get_transform<double>();

        // (1,0,0) after rotation about Z by 180° -> (-1,0,0), then + (1,1,0) = (0,1,0)
        CHECK(f({1,0,0}) == Vector3<double>(0,1,0));

        // (0,1,0) after rotation -> (0,-1,0), then translate = (1,0,0)
        CHECK(f({0,1,0}) == Vector3<double>(1,0,0));

        // (1,1,1) after rotation -> (-1,-1,1), then translate = (0,0,1)
        CHECK(f({1,1,1}) == Vector3<double>(0,0,1));
    }
}
