#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/parameters/RotationsOnly.h>

using namespace rigidbody;

TEST_CASE("RotationsOnly::next") {
    int iterations = 10;
    double length_start = GENERATE(1, 2, 3);
    double rad_start = GENERATE(1, 2, 3);
    rigidbody::parameter::RotationsOnly ro(iterations, length_start, rad_start);

    for (int i = 0; i < iterations; i++) {
        auto p = ro.next();
        REQUIRE(p.dr.x() == 0);
        REQUIRE(p.dr.x() == 0);
        REQUIRE(p.dr.y() == 0);
        REQUIRE(p.dr.y() == 0);
        REQUIRE(p.dr.z() == 0);
        REQUIRE(p.dr.z() == 0);
        REQUIRE(p.alpha >= -rad_start);
        REQUIRE(p.alpha <= rad_start);
        REQUIRE(p.beta >= -rad_start);
        REQUIRE(p.beta <= rad_start);
        REQUIRE(p.gamma >= -rad_start);
        REQUIRE(p.gamma <= rad_start);
    }
}
