#include <catch2/catch_test_macros.hpp>

#include <rigidbody/parameters/Parameter.h>

using namespace rigidbody;

TEST_CASE("Parameters::Parameter") {
    SECTION("default") {
        rigidbody::Parameter p;
        CHECK(p.dx == Vector3<double>(0, 0, 0));
        CHECK(p.alpha == 0);
        CHECK(p.beta == 0);
        CHECK(p.gamma == 0);
    }

    SECTION("Vector3<double>&, double, double, double") {
        Vector3<double> dx(1, 2, 3);
        double alpha = 4, beta = 5, gamma = 6;
        rigidbody::Parameter p(dx, alpha, beta, gamma);
        CHECK(p.dx == dx);
        CHECK(p.alpha == alpha);
        CHECK(p.beta == beta);
        CHECK(p.gamma == gamma);
    }
}