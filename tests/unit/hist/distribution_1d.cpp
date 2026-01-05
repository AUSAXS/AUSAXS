#include <catch2/catch_test_macros.hpp>

#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>

using namespace ausaxs;

TEST_CASE("Distribution1D::Distribution1D") {
    SECTION("default constructor") {
        hist::Distribution1D dist;
        CHECK(dist.size() == 0);
    }

    SECTION("size constructor") {
        hist::Distribution1D dist(10);
        CHECK(dist.size() == 10);
    }

    SECTION("vector constructor") {
        std::vector<double> data{1, 2, 3, 4, 5};
        hist::Distribution1D dist(data);
        CHECK(dist.size() == 5);
        CHECK(dist.as_vector() == data);
    }

    SECTION("from WeightedDistribution1D") {
        hist::WeightedDistribution1D wdist(10);
        wdist.set_content(0, 1);
        wdist.set_content(1, 2);
        wdist.set_content(2, 3);

        hist::Distribution1D dist(wdist);
        CHECK(dist.size() == 10);
        CHECK(dist.get_content(0) == 1);
        CHECK(dist.get_content(1) == 2);
        CHECK(dist.get_content(2) == 3);
    }
}

TEST_CASE("Distribution1D::get_content") {
    SECTION("const reference") {
        std::vector<double> data{1, 2, 3, 4, 5};
        hist::Distribution1D dist(data);
        CHECK(dist.get_content() == data);
    }

    SECTION("indexed access") {
        std::vector<double> data{1, 2, 3, 4, 5};
        hist::Distribution1D dist(data);
        CHECK(dist.get_content(0) == 1);
        CHECK(dist.get_content(2) == 3);
        CHECK(dist.get_content(4) == 5);
    }
}

TEST_CASE("Distribution1D::set_content") {
    hist::Distribution1D dist(5);
    dist.set_content(0, 10);
    dist.set_content(1, 20);
    dist.set_content(2, 30);
    CHECK(dist.get_content(0) == 10);
    CHECK(dist.get_content(1) == 20);
    CHECK(dist.get_content(2) == 30);
}

TEST_CASE("Distribution1D::add_index") {
    hist::Distribution1D dist(10);
    dist.add_index(0, 1);
    dist.add_index(1, 2);
    dist.add_index(2, 3);
    
    CHECK(dist.get_content(0) == 1);
    CHECK(dist.get_content(1) == 2);
    CHECK(dist.get_content(2) == 3);
}

TEST_CASE("Distribution1D::increment_index") {
    SECTION("basic increment") {
        hist::Distribution1D dist(10);
        dist.increment_index(0);
        dist.increment_index(1);
        dist.increment_index(1);
        dist.increment_index(2);
        dist.increment_index(2);
        dist.increment_index(2);
        
        CHECK(dist.get_content(0) == 1);
        CHECK(dist.get_content(1) == 2);
        CHECK(dist.get_content(2) == 3);
    }

    SECTION("template parameter increment") {
        hist::Distribution1D dist(10);
        dist.increment_index<2>(0);
        dist.increment_index<3>(1);
        dist.increment_index<5>(2);
        
        CHECK(dist.get_content(0) == 2);
        CHECK(dist.get_content(1) == 3);
        CHECK(dist.get_content(2) == 5);
    }
}

TEST_CASE("Distribution1D::increment_linear_index") {
    SECTION("basic increment") {
        hist::Distribution1D dist(10);
        dist.increment_linear_index(0);
        dist.increment_linear_index(1);
        dist.increment_linear_index(1);
        dist.increment_linear_index(2);
        dist.increment_linear_index(2);
        dist.increment_linear_index(2);
        
        CHECK(dist.get_content(0) == 1);
        CHECK(dist.get_content(1) == 2);
        CHECK(dist.get_content(2) == 3);
    }

    SECTION("template parameter increment") {
        hist::Distribution1D dist(10);
        dist.increment_linear_index<2>(0);
        dist.increment_linear_index<3>(1);
        dist.increment_linear_index<5>(2);
        
        CHECK(dist.get_content(0) == 2);
        CHECK(dist.get_content(1) == 3);
        CHECK(dist.get_content(2) == 5);
    }

    SECTION("linear_index equals index for 1D") {
        hist::Distribution1D dist(10);
        dist.increment_index(5);
        hist::Distribution1D dist2(10);
        dist2.increment_linear_index(5);
        
        CHECK(dist.get_content() == dist2.get_content());
    }
}

TEST_CASE("Distribution1D::clear") {
    hist::Distribution1D dist(5);
    dist.set_content(2, 10);
    CHECK(dist.get_content(2) == 10);
    
    dist.clear(2);
    CHECK(dist.get_content(2) == 0);
}

TEST_CASE("Distribution1D::as_vector") {
    std::vector<double> data{1, 2, 3, 4, 5};
    hist::Distribution1D dist(data);
    CHECK(dist.as_vector() == data);
}

TEST_CASE("Distribution1D::operator+=") {
    std::vector<double> data1{1, 2, 3, 4, 5};
    std::vector<double> data2{5, 4, 3, 2, 1};
    hist::Distribution1D dist1(data1);
    hist::Distribution1D dist2(data2);
    
    dist1 += dist2;
    CHECK(dist1.get_content() == std::vector<double>{6, 6, 6, 6, 6});
}

TEST_CASE("Distribution1D::operator-=") {
    std::vector<double> data1{5, 4, 3, 2, 1};
    std::vector<double> data2{1, 2, 3, 2, 1};
    hist::Distribution1D dist1(data1);
    hist::Distribution1D dist2(data2);
    
    dist1 -= dist2;
    CHECK(dist1.get_content() == std::vector<double>{4, 2, 0, 0, 0});
}

TEST_CASE("Distribution1D::operator*") {
    std::vector<double> data{1, 2, 3, 4, 5};
    hist::Distribution1D dist(data);
    
    auto result = 2.0 * dist;
    CHECK(result.get_content() == std::vector<double>{2, 4, 6, 8, 10});
}
