#include <catch2/catch_test_macros.hpp>

#include <utility/Container2D.h>

TEST_CASE("Container2D::Container2D") {
    SECTION("default") {
        Container2D<int> container;
        CHECK(container.N == 0);
        CHECK(container.M == 0);
    }

    SECTION("unsigned int, unsigned int, unsigned int") {
        Container2D<int> container(2, 3);
        CHECK(container.N == 2);
        CHECK(container.M == 3);
    }

    SECTION("unsigned int, unsigned int, unsigned int, const T&") {
        Container2D<int> container(2, 3, 5);
        CHECK(container.N == 2);
        CHECK(container.M == 3);
        for (unsigned int i = 0; i < container.N; ++i) {
            for (unsigned int j = 0; j < container.M; ++j) {
                CHECK(container(i, j) == 5);
            }
        }
    }
}

TEST_CASE("Container3D::iterators") {
    SECTION("default") {
        Container2D<int> container1(2, 3, 5);
        Container2D<int> container2(2, 3);
        std::copy(container1.begin(), container1.end(), container2.begin());
        for (unsigned int i = 0; i < container2.N; ++i) {
            for (unsigned int j = 0; j < container2.M; ++j) {
                CHECK(container2(i, j) == 5);
            }
        }
    }

    SECTION("index iterators") {
        SECTION("simple") {
            Container2D<int> container(2, 3, 5);
            std::vector<int> dest(3, 0);
            std::copy(container.begin(0), container.end(0), dest.begin());
            for (unsigned int i = 0; i < dest.size(); ++i) {
                CHECK(dest[i] == 5);
            }
        }

        SECTION("different vals") {
            Container2D<int> container(2, 3);
            int c = 0;
            for (unsigned int i = 0; i < container.N; ++i) {
                for (unsigned int j = 0; j < container.M; ++j) {
                    container(i, j) = c++;
                }
            }
            std::vector<int> dest(3, 0);
            std::copy(container.begin(0), container.end(0), dest.begin());
            CHECK(dest == std::vector<int>({0, 1, 2}));

            std::copy(container.begin(1), container.end(1), dest.begin());
            CHECK(dest == std::vector<int>({3, 4, 5}));
        }
    }
}