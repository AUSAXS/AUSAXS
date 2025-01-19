#include <catch2/catch_test_macros.hpp>

#include <rigidbody/parameters/Parameter.h>

using namespace ausaxs;
using namespace rigidbody;

TEST_CASE("Parameters::Parameter") {
    SECTION("default") {
        rigidbody::parameter::Parameter p;
        CHECK(p.translation == Vector3(0, 0, 0));
        CHECK(p.rotation == Vector3{0, 0, 0});
    }

    SECTION("Vector3<double>&, double, double, double") {
        Vector3<double> dx(1, 2, 3);
        Vector3<double> dr(4, 5, 6);
        rigidbody::parameter::Parameter p(dx, dr);
        CHECK(p.translation == dx);
        CHECK(p.rotation == dr);
    }
}