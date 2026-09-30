#include <catch2/catch_test_macros.hpp>

#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution3D.h>

#include <algorithm>

using namespace ausaxs;

TEST_CASE("Distribution3D::Distribution3D") {
    SECTION("default constructor") {
        hist::Distribution3D dist;
        CHECK(dist.size_x() == 0);
        CHECK(dist.size_y() == 0);
        CHECK(dist.size_z() == 0);
    }

    SECTION("size constructor") {
        hist::Distribution3D dist(5, 10, 15);
        CHECK(dist.size_x() == 5);
        CHECK(dist.size_y() == 10);
        CHECK(dist.size_z() == 15);
    }

    SECTION("from WeightedDistribution3D") {
        hist::WeightedDistribution3D wdist(5, 5, 10);
        wdist.index(0, 0, 0).value = 1;
        wdist.index(1, 1, 1).value = 2;
        wdist.index(2, 2, 2).value = 3;

        hist::Distribution3D dist(wdist);
        CHECK(dist.size_x() == 5);
        CHECK(dist.size_y() == 5);
        CHECK(dist.size_z() == 10);
        CHECK(dist.index(0, 0, 0) == 1);
        CHECK(dist.index(1, 1, 1) == 2);
        CHECK(dist.index(2, 2, 2) == 3);
    }
}

TEST_CASE("Distribution3D::add_index") {
    hist::Distribution3D dist(5, 5, 10);
    dist.add_index(0, 0, 0, 1);
    dist.add_index(1, 1, 1, 2);
    dist.add_index(2, 2, 2, 3);
    
    CHECK(dist.index(0, 0, 0) == 1);
    CHECK(dist.index(1, 1, 1) == 2);
    CHECK(dist.index(2, 2, 2) == 3);
}

TEST_CASE("Distribution3D::increment_linear_index") {
    SECTION("two parameters - ij, k access") {
        hist::Distribution3D dist(2, 2, 3);
        // increment_linear_index(ij, k) where ij is combined form factor index
        dist.increment_linear_index(0, 0);
        dist.increment_linear_index(0, 1);
        dist.increment_linear_index(1, 0);
        
        CHECK(dist.linear_index(0, 0) == 1);
        CHECK(dist.linear_index(0, 1) == 1);
        CHECK(dist.linear_index(1, 0) == 1);
    }

    SECTION("template parameter increment") {
        hist::Distribution3D dist(3, 3, 3);
        dist.increment_linear_index<2>(0, 0);
        dist.increment_linear_index<3>(0, 1);
        dist.increment_linear_index<5>(0, 2);
        
        CHECK(dist.linear_index(0, 0) == 2);
        CHECK(dist.linear_index(0, 1) == 3);
        CHECK(dist.linear_index(0, 2) == 5);
    }
}

TEST_CASE("Distribution3D: triangular") {
    constexpr auto Tri = hist::Shape::Triangular;
    auto fill = [] (auto& d) {
        for (int i = 0; i < d.size_x(); ++i) {for (int j = i; j < d.size_y(); ++j) {for (int k = 0; k < d.size_z(); ++k) {
            d.add_index(i, j, k, 0.25 + 0.5*k, 100*i + 10*j + k + 1); // distance, weight
        }}}
    };

    hist::WeightedDistribution3D<Tri> tri(4, 4, 3);
    hist::WeightedDistribution3D<> square(4, 4, 3);
    fill(tri); fill(square);
    CHECK(std::distance(tri.begin(), tri.end()) == 10*3);
    CHECK(&tri.index(3, 1, 2) == &tri.index(1, 3, 2));
    for (int i = 0; i < 4; ++i) {for (int j = 0; j < 4; ++j) {for (int k = 0; k < 3; ++k) {
        CHECK(tri.index(i, j, k).value == square.index(std::min(i, j), std::max(i, j), k).value);
    }}}

    // a square distribution filled on one half only is what the managers produced before, so the two must agree
    CHECK(tri.get_weights() == square.get_weights());
    hist::Distribution3D<Tri> plain(tri);
    plain.resize(2);
    CHECK(std::distance(plain.begin(), plain.end()) == 10*2);
    CHECK(plain.index(2, 1, 1) == 100*1 + 10*2 + 1 + 1);
    CHECK(hist::WeightedDistribution3D<Tri>(plain).index(1, 2, 1).value == plain.index(2, 1, 1));
}
