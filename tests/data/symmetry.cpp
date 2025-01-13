#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Symmetry.h>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("Symmetry::get_transform") {
    SECTION("no transformation") {
        // Default constructor: no translation, no rotation
        data::detail::Symmetry s;
        auto f = s.get_transform<double>({0, 0, 0});
        CHECK(f({1, 2, 3}) == Vector3<double>(1, 2, 3));
    }

    SECTION("translation only") {
        // Translate by (1,2,3)
        data::detail::Symmetry s({1, 2, 3});
        auto f = s.get_transform<double>({0, 0, 0});
        CHECK(f({1, 2, 3}) == Vector3<double>(2, 4, 6));
        CHECK(f({0, 0, 0}) == Vector3<double>(1, 2, 3));
        CHECK(f({-1, -1, -1}) == Vector3<double>(0, 1, 2));
    }

    SECTION("internal rotation about X by +90 deg") {
        // A rotation about X by +90° (π/2): 
        // (x, y, z) -> (x, z, -y)
        data::detail::Symmetry s({0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {std::numbers::pi/2, 0, 0}, 1);
        auto f = s.get_transform<double>({0, 0, 0});
        
        CHECK(f({1, 0, 0}) == Vector3<double>(1,  0, 0));   // x-axis unchanged
        CHECK(f({0, 1, 0}) == Vector3<double>(0,  0, 1));   // y -> z
        CHECK(f({0, 0, 1}) == Vector3<double>(0, -1, 0));   // z -> -y
    }

    SECTION("external rotation about Y by +90 deg") {
        // Rotate about Y by +90° (π/2):
        // (x, y, z) -> (z, y, -x)
        data::detail::Symmetry s({0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/2, 0}, {0, 0, 0}, 1);
        auto f = s.get_transform<double>({0, 0, 0});
        
        CHECK(f({1, 0, 0}) == Vector3<double>(0, 0, -1));  // x -> -z
        CHECK(f({0, 1, 0}) == Vector3<double>(0, 1, 0));   // y unchanged
        CHECK(f({0, 0, 1}) == Vector3<double>(1, 0, 0));   // z -> x
    }

    SECTION("translation + internal rotation about X by +90 deg") {
        // Combine a translation (1,2,3) and internal rotation about X by +90°
        // Rotation about X by +90°: (x, y, z) -> (x, z, -y)
        // Then translate by (1,2,3).
        data::detail::Symmetry s({1, 2, 3}, {0, 0, 0}, {0, 0, 0}, {std::numbers::pi/2, 0, 0}, 1);
        auto f = s.get_transform<double>({0, 0, 0});
        
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

    SECTION("external rotation about Z by 180 deg + translation") {
        // Rotate about Z by 180°: (x, y, z) -> (-x, -y, z)
        // Then translate by (1,1,0).
        data::detail::Symmetry s({1, 1, 0}, {0, 0, 0}, {0, 0, std::numbers::pi}, {0,0,0}, 1);
        auto f = s.get_transform<double>({0, 0, 0});

        // (1,0,0) after rotation about Z by 180° -> (-1,0,0), then + (1,1,0) = (0,1,0)
        CHECK(f({1,0,0}) == Vector3<double>(0,1,0));

        // (0,1,0) after rotation -> (0,-1,0), then translate = (1,0,0)
        CHECK(f({0,1,0}) == Vector3<double>(1,0,0));

        // (1,1,1) after rotation -> (-1,-1,1), then translate = (0,0,1)
        CHECK(f({1,1,1}) == Vector3<double>(0,0,1));
    }
}
