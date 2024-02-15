#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/parameters/SimpleParameterGeneration.h>

using namespace rigidbody;

TEST_CASE("SimpleParameterGeneration::next") {
    int iterations = 10;
    double length_start = GENERATE(1, 2, 3);
    double rad_start = GENERATE(1, 2, 3);
    rigidbody::parameter::SimpleParameterGeneration spg(iterations, length_start, rad_start);

    for (int i = 0; i < iterations; i++) {
        auto p = spg.next();
        REQUIRE(p.dr.x() >= -length_start);
        REQUIRE(p.dr.x() <= length_start);
        REQUIRE(p.dr.y() >= -length_start);
        REQUIRE(p.dr.y() <= length_start);
        REQUIRE(p.dr.z() >= -length_start);
        REQUIRE(p.dr.z() <= length_start);
        REQUIRE(p.alpha >= -rad_start);
        REQUIRE(p.alpha <= rad_start);
        REQUIRE(p.beta >= -rad_start);
        REQUIRE(p.beta <= rad_start);
        REQUIRE(p.gamma >= -rad_start);
        REQUIRE(p.gamma <= rad_start);
    }
}
