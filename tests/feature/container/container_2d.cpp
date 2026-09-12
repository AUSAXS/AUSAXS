#include <algorithm>
#include <catch2/catch_test_macros.hpp>

#include <container/Container2D.h>


using namespace ausaxs;
using namespace container;

TEST_CASE("Container2D::Container2D") {
    SECTION("default") {
        Container2D<int> container;
        CHECK(container.size_x() == 0);
        CHECK(container.size_y() == 0);
    }

    SECTION("int, int, int") {
        Container2D<int> container(2, 3);
        CHECK(container.size_x() == 2);
        CHECK(container.size_y() == 3);
    }

    SECTION("int, int, int, const T&") {
        Container2D<int> container(2, 3, 5);
        CHECK(container.size_x() == 2);
        CHECK(container.size_y() == 3);
        for (int i = 0; i < container.size_x(); ++i) {
            for (int j = 0; j < container.size_y(); ++j) {
                CHECK(container(i, j) == 5);
            }
        }
    }
}

TEST_CASE("Container3D::iterators") {
    SECTION("default") {
        Container2D<int> container1(2, 3, 5);
        Container2D<int> container2(2, 3);
        std::ranges::copy(container1, container2.begin());
        for (int i = 0; i < container2.size_x(); ++i) {
            for (int j = 0; j < container2.size_y(); ++j) {
                CHECK(container2(i, j) == 5);
            }
        }
    }

    SECTION("index iterators") {
        SECTION("simple") {
            Container2D<int> container(2, 3, 5);
            std::vector<int> dest(3, 0);
            std::copy(container.begin(0), container.end(0), dest.begin());
            for (int i : dest) {
                CHECK(i == 5);
            }
        }

        SECTION("different vals") {
            Container2D<int> container(2, 3);
            int c = 0;
            for (int i = 0; i < container.size_x(); ++i) {
                for (int j = 0; j < container.size_y(); ++j) {
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

TEST_CASE("Container2D::resize") {
    Container2D<double> container(3, 4);
    {
        int i = 0;
        std::ranges::transform(container, container.begin(), [&i](double) {return i++;});
    }

    SECTION("larger") {
        container.resize(5);
        CHECK(container.size_x() == 3);
        CHECK(container.size_y() == 5);

        int c = 0;
        for (int i = 0; i < container.size_x(); ++i) {
            for (int j = 0; j < 4; ++j) {
                CHECK(container(i, j) == c++);
            }
            CHECK(container(i, 4) == 0);
        }
    }

    SECTION("smaller") {
        container.resize(3);
        CHECK(container.size_x() == 3);
        CHECK(container.size_y() == 3);

        int c = 0;
        for (int i = 0; i < container.size_x(); ++i) {
            for (int j = 0; j < container.size_y(); ++j) {
                CHECK(container(i, j) == c++);
            }
            c++;
        }
    }
}