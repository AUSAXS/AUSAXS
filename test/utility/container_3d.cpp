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
    SECTION("default") {
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

    SECTION("index iterators") {
        SECTION("simple") {
            Container3D<int> container(2, 3, 4, 5);
            std::vector<int> dest(4, 0);
            std::copy(container.begin(0, 0), container.end(0, 0), dest.begin());
            for (unsigned int i = 0; i < dest.size(); ++i) {
                CHECK(dest[i] == 5);
            }
        }

        SECTION("different vals") {
            Container3D<int> container(2, 3, 4);
            int c = 0;
            for (unsigned int i = 0; i < container.N; ++i) {
                for (unsigned int j = 0; j < container.M; ++j) {
                    for (unsigned int k = 0; k < container.L; ++k) {
                        container(i, j, k) = c++;
                    }
                }
            }
            std::vector<int> dest(4, 0);

            std::copy(container.begin(0, 0), container.end(0, 0), dest.begin());
            CHECK(dest == std::vector<int>({0, 1, 2, 3}));

            std::copy(container.begin(0, 1), container.end(0, 1), dest.begin());
            CHECK(dest == std::vector<int>({4, 5, 6, 7}));

            std::copy(container.begin(0, 2), container.end(0, 2), dest.begin());
            CHECK(dest == std::vector<int>({8, 9, 10, 11}));

            std::copy(container.begin(1, 0), container.end(1, 0), dest.begin());
            CHECK(dest == std::vector<int>({12, 13, 14, 15}));

            std::copy(container.begin(1, 1), container.end(1, 1), dest.begin());
            CHECK(dest == std::vector<int>({16, 17, 18, 19}));

            std::copy(container.begin(1, 2), container.end(1, 2), dest.begin());
            CHECK(dest == std::vector<int>({20, 21, 22, 23}));
        }
    }
}