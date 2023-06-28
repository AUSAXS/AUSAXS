#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/Image.h>

TEST_CASE("Image::Image") {
    SECTION("std::shared_ptr<ccp4::Header>, unsigned int") {}
    SECTION("Matrix<float>&") {}
    SECTION("Matrix<float>&, std::shared_ptr<ccp4::Header>, unsigned int") {}
    CHECK(false);
}

TEST_CASE("Image::as_hist") {}
TEST_CASE("Image::generate_atoms") {}
TEST_CASE("Image::count_voxels") {}
TEST_CASE("Image::set_z") {}
TEST_CASE("Image::get_z") {}
TEST_CASE("Image::mean") {}
TEST_CASE("Image::limits") {}
TEST_CASE("Image::get_bounds") {}
TEST_CASE("Image::set_header") {}
TEST_CASE("Image::setup_bounds") {}
TEST_CASE("Image::index") {}
TEST_CASE("Image::squared_sum") {}
TEST_CASE("Image::operator==") {}
TEST_CASE("Image::to_string") {}