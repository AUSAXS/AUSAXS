#include <catch2/catch_test_macros.hpp>

#include <hist/detail/CompactCoordinatesData.h>
#include <math/Vector3.h>

TEST_CASE("hist::detail::CompactCoordinatesData") {
    SECTION("default") {
        hist::detail::CompactCoordinatesData data;
        CHECK(data.x == 0);
        CHECK(data.y == 0);
        CHECK(data.z == 0);
        CHECK(data.w == 0);
    }

    SECTION("Vector3<double>, float") {
        hist::detail::CompactCoordinatesData data(Vector3<double>(1, 2, 3), 4);
        CHECK(data.x == 1);
        CHECK(data.y == 2);
        CHECK(data.z == 3);
        CHECK(data.w == 4);
    }

    SECTION("evaluate") {
        hist::detail::CompactCoordinatesData data1({1, 1, 1}, 2);
        hist::detail::CompactCoordinatesData data2({2, 1, 1}, 4);
        auto result = data1.evaluate(data2);
        CHECK(result.distance == 1);
        CHECK(result.weight == 8);
    }
}