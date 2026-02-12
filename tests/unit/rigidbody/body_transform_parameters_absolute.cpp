#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>
#include <math/MatrixUtils.h>
#include <numbers>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

TEST_CASE("BodyTransformParametersAbsolute::transform with pivot and rotation") {
    SECTION("rotation around origin") {
        BodyTransformParametersAbsolute params({0, 0, 0}, {0, 0, 0});
        Vector3<double> pivot(0, 0, 0);
        Vector3<double> euler_angles(0, 0, std::numbers::pi/2);
        auto rotation_matrix = matrix::rotation_matrix(euler_angles);
        
        params.transform(pivot, rotation_matrix);
        
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(params.rotation.z(), Catch::Matchers::WithinAbs(std::numbers::pi/2, 1e-6));
    }

    SECTION("rotation around non-origin pivot") {
        BodyTransformParametersAbsolute params({5, 0, 0}, {0, 0, 0});
        Vector3<double> pivot(0, 0, 0);
        Vector3<double> euler_angles(0, 0, std::numbers::pi/2);
        auto rotation_matrix = matrix::rotation_matrix(euler_angles);
        
        params.transform(pivot, rotation_matrix);
        
        // After 90 degree rotation around z-axis, (5,0,0) becomes (0,5,0)
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(5.0, 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(0.0, 1e-10));
    }

    SECTION("rotation with offset pivot") {
        BodyTransformParametersAbsolute params({3, 0, 0}, {0, 0, 0});
        Vector3<double> pivot(1, 0, 0);
        Vector3<double> euler_angles(0, 0, std::numbers::pi/2);
        auto rotation_matrix = matrix::rotation_matrix(euler_angles);
        
        params.transform(pivot, rotation_matrix);
        
        // Translation is at (3,0,0), pivot is at (1,0,0)
        // Relative position is (2,0,0), after rotation becomes (0,2,0)
        // Final position is pivot + rotated = (1,0,0) + (0,2,0) = (1,2,0)
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(0.0, 1e-10));
    }

    SECTION("rotation accumulates with existing rotation") {
        BodyTransformParametersAbsolute params({0, 0, 0}, {0, 0, std::numbers::pi/4});
        Vector3<double> pivot(0, 0, 0);
        Vector3<double> euler_angles(0, 0, std::numbers::pi/4);
        auto rotation_matrix = matrix::rotation_matrix(euler_angles);
        
        params.transform(pivot, rotation_matrix);
        
        // Two 45-degree rotations should give 90 degrees
        REQUIRE_THAT(params.rotation.z(), Catch::Matchers::WithinAbs(std::numbers::pi/2, 1e-6));
    }

    SECTION("180 degree rotation") {
        BodyTransformParametersAbsolute params({1, 0, 0}, {0, 0, 0});
        Vector3<double> pivot(0, 0, 0);
        Vector3<double> euler_angles(0, 0, std::numbers::pi);
        auto rotation_matrix = matrix::rotation_matrix(euler_angles);
        
        params.transform(pivot, rotation_matrix);
        
        // After 180 degree rotation, (1,0,0) becomes (-1,0,0)
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(-1.0, 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("BodyTransformParametersAbsolute::transform with translation") {
    SECTION("simple translation") {
        BodyTransformParametersAbsolute params({0, 0, 0}, {0, 0, 0});
        Vector3<double> delta(1, 2, 3);
        
        params.transform(delta);
        
        CHECK(params.translation.x() == 1);
        CHECK(params.translation.y() == 2);
        CHECK(params.translation.z() == 3);
        CHECK(params.rotation.x() == 0);
        CHECK(params.rotation.y() == 0);
        CHECK(params.rotation.z() == 0);
    }

    SECTION("accumulated translations") {
        BodyTransformParametersAbsolute params({1, 2, 3}, {0, 0, 0});
        Vector3<double> delta(4, 5, 6);
        
        params.transform(delta);
        
        CHECK(params.translation.x() == 5);
        CHECK(params.translation.y() == 7);
        CHECK(params.translation.z() == 9);
    }

    SECTION("negative translation") {
        BodyTransformParametersAbsolute params({10, 10, 10}, {0, 0, 0});
        Vector3<double> delta(-5, -5, -5);
        
        params.transform(delta);
        
        CHECK(params.translation.x() == 5);
        CHECK(params.translation.y() == 5);
        CHECK(params.translation.z() == 5);
    }
}

TEST_CASE("BodyTransformParametersAbsolute::transform combined") {
    SECTION("rotation then translation") {
        BodyTransformParametersAbsolute params({1, 0, 0}, {0, 0, 0});
        Vector3<double> pivot(0, 0, 0);
        Vector3<double> euler_angles(0, 0, std::numbers::pi/2);
        auto rotation_matrix = matrix::rotation_matrix(euler_angles);
        Vector3<double> translation(1, 1, 1);
        
        params.transform(pivot, rotation_matrix, translation);
        
        // First rotation: (1,0,0) -> (0,1,0)
        // Then translation: (0,1,0) + (1,1,1) = (1,2,1)
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(1.0, 1e-10));
    }

    SECTION("multiple transformations preserve order") {
        BodyTransformParametersAbsolute params({0, 0, 0}, {0, 0, 0});
        Vector3<double> pivot(0, 0, 0);
        
        // Apply first rotation
        Vector3<double> euler1(0, 0, std::numbers::pi/4);
        auto rot1 = matrix::rotation_matrix(euler1);
        params.transform(pivot, rot1);
        
        // Apply translation
        params.transform(Vector3<double>(1, 0, 0));
        
        // Apply second rotation
        Vector3<double> euler2(0, 0, std::numbers::pi/4);
        auto rot2 = matrix::rotation_matrix(euler2);
        params.transform(pivot, rot2);
        
        // Combined rotations should be pi/2
        REQUIRE_THAT(params.rotation.z(), Catch::Matchers::WithinAbs(std::numbers::pi/2, 1e-6));
    }
}

TEST_CASE("BodyTransformParametersAbsolute::transform identity operations") {
    SECTION("identity rotation does not change translation") {
        BodyTransformParametersAbsolute params({1, 2, 3}, {0.1, 0.2, 0.3});
        Vector3<double> pivot(0, 0, 0);
        Vector3<double> euler_zero(0, 0, 0);
        auto identity = matrix::rotation_matrix(euler_zero);
        
        auto original_translation = params.translation;
        params.transform(pivot, identity);
        
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(original_translation.x(), 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(original_translation.y(), 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(original_translation.z(), 1e-10));
    }

    SECTION("zero translation does not change position") {
        BodyTransformParametersAbsolute params({1, 2, 3}, {0.1, 0.2, 0.3});
        Vector3<double> zero(0, 0, 0);
        
        auto original = params.translation;
        params.transform(zero);
        
        CHECK(params.translation.x() == original.x());
        CHECK(params.translation.y() == original.y());
        CHECK(params.translation.z() == original.z());
    }
}

TEST_CASE("BodyTransformParametersAbsolute::transform edge cases") {
    SECTION("large rotation angles") {
        BodyTransformParametersAbsolute params({1, 0, 0}, {0, 0, 0});
        Vector3<double> pivot(0, 0, 0);
        
        // Full rotation (2*pi)
        Vector3<double> euler_full(0, 0, 2*std::numbers::pi);
        auto rotation = matrix::rotation_matrix(euler_full);
        params.transform(pivot, rotation);
        
        // Should be back to original position
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(0.0, 1e-10));
    }

    SECTION("3D rotation around all axes") {
        BodyTransformParametersAbsolute params({1, 1, 1}, {0, 0, 0});
        Vector3<double> pivot(0, 0, 0);
        
        // Rotation around all three axes
        Vector3<double> euler_angles(0.1, 0.2, 0.3);
        auto rotation = matrix::rotation_matrix(euler_angles);
        params.transform(pivot, rotation);
        
        // Just verify no crashes and rotation is applied
        REQUIRE(params.rotation.x() != 0);
        REQUIRE(params.rotation.y() != 0);
        REQUIRE(params.rotation.z() != 0);
    }
}
