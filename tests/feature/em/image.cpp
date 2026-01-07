#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/Image.h>
#include <em/ImageStack.h>
#include <em/ObjectBounds3D.h>
#include <numeric>
#include <utility/Limit.h>

using namespace ausaxs;

// Unit tests for Image have been moved to tests/unit/em/image.cpp
// This file now contains only feature/integration tests

// TODO: fix the remaining tests
TEST_CASE("Image::Image") {
    SECTION("std::shared_ptr<ccp4::Header>, unsigned int") {}
    SECTION("Matrix<float>&") {}
    SECTION("Matrix<float>&, std::shared_ptr<ccp4::Header>, unsigned int") {}
}

TEST_CASE("Image::as_hist") {}
TEST_CASE("Image::generate_atoms") {}
TEST_CASE("Image::set_header") {}

TEST_CASE("Image::setup_bounds") {
    SECTION("correct_bounds_imagestack") {
        Matrix data = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
        Matrix data2 = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
        em::ImageStack images({data, data2});
        
        em::ObjectBounds3D bounds = images.minimum_volume(1);
        CHECK(bounds.total_volume() == 2*6*6);
        CHECK(bounds.bounded_volume() == ((4 + 3 + 3 + 4 + 3 + 5) + (5 + 4 + 4 + 5 + 4 + 6)));

        bounds = images.minimum_volume(2);
        CHECK(bounds.bounded_volume() == ((2 + 3 + 2 + 3 + 2 + 3) + (3 + 4 + 1 + 4 + 1 + 4)));
    }
}