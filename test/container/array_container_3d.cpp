#include <catch2/catch_test_macros.hpp>

#include <container/ArrayContainer3D.h>

container::ArrayContainer3D<double, 1, 2, 3> dummy;
TEST_CASE("ArrayContainer3D::ArrayContainer3D") {
    SECTION("default") {
        CHECK(dummy.N == 1);
        CHECK(dummy.M == 2);
        CHECK(dummy.L == 3);
    }
}