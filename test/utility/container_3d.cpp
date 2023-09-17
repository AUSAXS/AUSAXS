#include <catch2/catch_test_macros.hpp>

#include <utility/Container3D.h>

TEST_CASE("Container3D::Container3D") {
    SECTION("default") {
        Container3D<int> container;
        CHECK(container.N == 0);
        CHECK(container.M == 0);
        CHECK(container.L == 0);
    }

    SECTION("unsigned int, unsigned int, unsigned int") {
        Container3D<int> container(2, 3, 4);
        CHECK(container.N == 2);
        CHECK(container.M == 3);
        CHECK(container.L == 4);
    }

    SECTION("unsigned int, unsigned int, unsigned int, const T&") {
        Container3D<int> container(2, 3, 4, 5);
        CHECK(container.N == 2);
        CHECK(container.M == 3);
        CHECK(container.L == 4);
        for (unsigned int i = 0; i < container.N; ++i) {
            for (unsigned int j = 0; j < container.M; ++j) {
                for (unsigned int k = 0; k < container.L; ++k) {
                    CHECK(container(i, j, k) == 5);
                }
            }
        }
    }
}

TEST_CASE("Container3D::iterators") {
    Container3D<int> container1(2, 3, 4, 5);
    Container3D<int> container2(2, 3, 4);
    std::copy(container1.begin(), container1.end(), container2.begin());
    for (unsigned int i = 0; i < container2.N; ++i) {
        for (unsigned int j = 0; j < container2.M; ++j) {
            for (unsigned int k = 0; k < container2.L; ++k) {
                CHECK(container2(i, j, k) == 5);
            }
        }
    }
}