#include <catch2/catch_test_macros.hpp>

#include <container/Container3D.h>

using namespace container;

TEST_CASE("Container3D::Container3D") {
    SECTION("default") {
        Container3D<int> container;
        CHECK(container.size_x() == 0);
        CHECK(container.size_y() == 0);
        CHECK(container.size_z() == 0);
    }

    SECTION("unsigned int, unsigned int, unsigned int") {
        Container3D<int> container(2, 3, 4);
        CHECK(container.size_x() == 2);
        CHECK(container.size_y() == 3);
        CHECK(container.size_z() == 4);
    }

    SECTION("unsigned int, unsigned int, unsigned int, const T&") {
        Container3D<int> container(2, 3, 4, 5);
        CHECK(container.size_x() == 2);
        CHECK(container.size_y() == 3);
        CHECK(container.size_z() == 4);
        for (unsigned int i = 0; i < container.size_x(); ++i) {
            for (unsigned int j = 0; j < container.size_y(); ++j) {
                for (unsigned int k = 0; k < container.size_z(); ++k) {
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
        for (unsigned int i = 0; i < container2.size_x(); ++i) {
            for (unsigned int j = 0; j < container2.size_y(); ++j) {
                for (unsigned int k = 0; k < container2.size_z(); ++k) {
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
            for (unsigned int i = 0; i < container.size_x(); ++i) {
                for (unsigned int j = 0; j < container.size_y(); ++j) {
                    for (unsigned int k = 0; k < container.size_z(); ++k) {
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

TEST_CASE("Container3D::resize") {
    Container3D<double> container(2, 3, 4);
    {
        unsigned int i = 0;
        std::transform(container.begin(), container.end(), container.begin(), [&i](double) {return i++;});
    }

    SECTION("larger") {
        container.resize(5);
        CHECK(container.size_x() == 2);
        CHECK(container.size_y() == 3);
        CHECK(container.size_z() == 5);

        unsigned int c = 0;
        for (unsigned int i = 0; i < container.size_x(); ++i) {
            for (unsigned int j = 0; j < container.size_y(); ++j) {
                for (unsigned int k = 0; k < 4; ++k) {
                    CHECK(container(i, j, k) == c++);
                }
                CHECK(container(i, j, 4) == 0);
            }
        }
    }

    SECTION("smaller") {
        container.resize(3);
        CHECK(container.size_x() == 2);
        CHECK(container.size_y() == 3);
        CHECK(container.size_z() == 3);

        unsigned int c = 0;
        for (unsigned int i = 0; i < container.size_x(); ++i) {
            for (unsigned int j = 0; j < container.size_y(); ++j) {
                for (unsigned int k = 0; k < container.size_z(); ++k) {
                    CHECK(container(i, j, k) == i*3*4 + j*4 + k);
                }
            }
            c--;
        }
    }
}