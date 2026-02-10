#include <catch2/catch_test_macros.hpp>

#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::parameter;

TEST_CASE("BodyTransformParametersRelative::BodyTransformParametersRelative") {
    SECTION("default") {
        BodyTransformParametersRelative params;
        CHECK_FALSE(params.translation.has_value());
        CHECK_FALSE(params.rotation.has_value());
        CHECK_FALSE(params.symmetry_pars.has_value());
    }

    SECTION("with translation and rotation") {
        Vector3<double> t = {1, 2, 3};
        Vector3<double> r = {0.1, 0.2, 0.3};
        BodyTransformParametersRelative params(t, r);
        REQUIRE(params.translation.has_value());
        REQUIRE(params.rotation.has_value());
        CHECK(params.translation.value() == t);
        CHECK(params.rotation.value() == r);
    }

    SECTION("equality") {
        BodyTransformParametersRelative a({1, 2, 3}, {0.1, 0.2, 0.3});
        BodyTransformParametersRelative b({1, 2, 3}, {0.1, 0.2, 0.3});
        BodyTransformParametersRelative c({1, 2, 4}, {0.1, 0.2, 0.3});
        CHECK(a == b);
        CHECK_FALSE(a == c);
    }
}

TEST_CASE("BodyTransformParametersAbsolute::BodyTransformParametersAbsolute") {
    SECTION("with translation and rotation") {
        Vector3<double> t = {1, 2, 3};
        Vector3<double> r = {0.1, 0.2, 0.3};
        BodyTransformParametersAbsolute params(t, r);
        CHECK(params.translation == t);
        CHECK(params.rotation == r);
    }
}
