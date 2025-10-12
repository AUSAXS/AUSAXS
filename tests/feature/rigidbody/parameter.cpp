#include <catch2/catch_test_macros.hpp>

#include <rigidbody/parameters/Parameter.h>

using namespace ausaxs;
using namespace rigidbody;

TEST_CASE("Parameters::Parameter") {
    SECTION("Vector3<double>, Vector3<double>") {
        Vector3<double> dx(1, 2, 3);
        Vector3<double> dr(4, 5, 6);
        rigidbody::parameter::Parameter p(dx, dr);
        CHECK(p.translation == dx);
        CHECK(p.rotation == dr);
    }
}