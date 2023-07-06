#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <vector>
#include <string>
#include <iostream>

#include <math/Matrix.h>
#include <math/Vector.h>
#include <math/Vector3.h>

// Slices are not meant to be used directly, but rather through the Matrix and Vector classes.
// Thus we will test them through these two. 
TEST_CASE("Slice::Slice") {
    Matrix A = {{1, 1, 2, 2}, {3, 3, 2, 2}, {5, 5, 4, 4}};

    SECTION("rows") {
        auto a = A.row(0);
        REQUIRE(a.size() == 4);
        CHECK(a == std::vector{1, 1, 2, 2});
        CHECK(a == A[0]);

        auto b = A.row(1);
        REQUIRE(b.size() == 4);
        CHECK(b == std::vector{3, 3, 2, 2});
        CHECK(b == A[1]);

        auto c = A.row(2);
        REQUIRE(c.size() == 4);
        CHECK(c == std::vector{5, 5, 4, 4});
        CHECK(c == A[2]);        
    }

    SECTION("cols") {
        auto a = A.col(0);
        REQUIRE(a.size() == 3);
        CHECK(a == std::vector{1, 3, 5});

        auto b = A.col(1);
        REQUIRE(b.size() == 3);
        CHECK(b == std::vector{1, 3, 5});

        auto c = A.col(2);
        REQUIRE(c.size() == 3);
        CHECK(c == std::vector{2, 2, 4});

        auto d = A.col(3);
        REQUIRE(d.size() == 3);
        CHECK(d == std::vector{2, 2, 4});
    }
}

TEST_CASE("Slice::operator=") {
    Matrix A = {{1, 1, 2, 2}, {3, 3, 2, 2}, {5, 5, 4, 4}};

    // row assignment
    A.row(1) = {9, 1, 2, 3};
    A.row(2) = {6, 3, 1, 2};
    REQUIRE(A == Matrix{{1, 1, 2, 2}, {9, 1, 2, 3}, {6, 3, 1, 2}});

    // column assignment
    A.col(1) = {2, 5, 1};
    A.col(3) = {7, 1, 3};
    REQUIRE(A == Matrix{{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}});

    // minus-assignment
    A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
    A.row(0) -= A.row(1);
    A.row(1) -= A.row(2);
    REQUIRE(A.row(0) == Vector<double>{-8, -3, 0, 6});
    REQUIRE(A.row(1) == Vector<double>{3, 4, 1, -2});

    A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
    A.col(0) -= A.col(1);
    A.col(1) -= A.col(2);
    REQUIRE(A.col(0) == Vector<double>{-1, 4, 5});
    REQUIRE(A.col(1) == Vector<double>{0, 3, 0});

    // plus-assignment
    A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
    A.row(0) += A.row(1);
    A.row(1) += A.row(2);
    REQUIRE(A.row(0) == Vector<double>{10, 7, 4, 8});
    REQUIRE(A.row(1) == Vector<double>{15, 6, 3, 4});

    A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
    A.col(0) += A.col(1);
    A.col(1) += A.col(2);
    REQUIRE(A.col(0) == Vector<double>{3, 14, 7});
    REQUIRE(A.col(1) == Vector<double>{4, 7, 2});
}

TEST_CASE("copy") {
    Matrix<double> A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
    auto a = A.row(1);

    auto b = a;
    REQUIRE(b == Vector<double>{9, 5, 2, 1});
    REQUIRE(b[0] == 9);
    REQUIRE(b[1] == 5);
    REQUIRE(b[2] == 2);
    REQUIRE(b[3] == 1);
}

TEST_CASE("vector_cast") {
    Matrix<double> A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
    Vector<double> a = A.col(2);
    REQUIRE((a.N == 3 && a.data.size() == 3));
    REQUIRE(a == Vector{2, 2, 1});
    // REQUIRE((A.col(2).operator Vector<double>().N == 3 && A.col(2).operator Vector<double>().data.size() == 3)); // chain cast
}

TEST_CASE("dot_product") {
    Matrix<double> A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};

    // dot with vector
    Vector<double> b = {2, 3, 1, 5};
    REQUIRE(A.row(0).dot(b) == (2+6+2+35));
    REQUIRE(A.row(2).dot(b) == (12+3+1+15));

    b = {1, 4, 2};
    REQUIRE(A.col(0).dot(b) == (1+36+12));
    REQUIRE(A.col(2).dot(b) == (2+8+2));

    // dot with other slice
    REQUIRE(A.col(0).dot(A.col(2)) == (2+18+6));
    REQUIRE(A.row(0).dot(A.row(1)) == (9+10+4+7));

}

TEST_CASE("norm") {
    Matrix<double> A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
    REQUIRE(A.col(0).norm() == sqrt(1+81+36));
    REQUIRE(A.row(0).norm() == sqrt(1+4+4+49));
}

TEST_CASE("iterators") {
    SECTION("basic") {
        Matrix<double> A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
        auto r = A.row(0);
        auto c = A.col(0);

        // basic iterator stuff
        REQUIRE(r.begin() == r.begin());
        REQUIRE(r.end() == r.end());
        REQUIRE(c.begin() == c.begin());
        REQUIRE(c.end() == c.end());

        REQUIRE(r.begin() != r.end());
        REQUIRE(c.begin() != c.end());

        // row iterator
        auto it = r.begin();
        CHECK(*it == 1);
        ++it;
        CHECK(*it == 2);
        ++it;
        CHECK(*it == 2);
        ++it;
        CHECK(*it == 7);

        for (auto i : r) {
            CHECK((i == 1 || i == 2 || i == 7));
        }

        // column iterator
        it = c.begin();
        CHECK(*it == 1);
        ++it;
        CHECK(*it == 9);
        ++it;
        CHECK(*it == 6);

        for (auto i : c) {
            CHECK((i == 1 || i == 9 || i == 6));
        }
    }

    SECTION("advanced") {
        Matrix<double> A = {{1, 2, 2, 7}, {9, 5, 2, 1}, {6, 1, 1, 3}};
        SECTION("single transform") {
            auto r = A.row(0);
            std::transform(r.begin(), r.end(), r.begin(), [](double x) { return x * 2; });
            REQUIRE(r == Vector{2, 4, 4, 14});
            REQUIRE(A.row(0) == Vector{2, 4, 4, 14});
        }

        SECTION("double transform") {
            auto r1 = A.row(0);
            auto r2 = A.row(1);
            std::transform(r1.begin(), r1.end(), r2.begin(), r1.begin(), std::plus<double>());
            REQUIRE(r1 == Vector{10, 7, 4, 8});
            REQUIRE(A.row(0) == Vector{10, 7, 4, 8});
        }
    }
}