#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/Image.h>
#include <em/ImageStack.h>
#include <em/ObjectBounds3D.h>
#include <utility/Limit.h>

TEST_CASE("Image::Image") {
    SECTION("std::shared_ptr<ccp4::Header>, unsigned int") {}
    SECTION("Matrix<float>&") {}
    SECTION("Matrix<float>&, std::shared_ptr<ccp4::Header>, unsigned int") {}
    CHECK(false);
}

TEST_CASE("Image::as_hist") {
    CHECK(false);
}

TEST_CASE("Image::generate_atoms") {
    CHECK(false);
}

TEST_CASE("Image::count_voxels") {
    CHECK(false);
}

TEST_CASE("Image::set_z") {
    CHECK(false);
}

TEST_CASE("Image::get_z") {
    CHECK(false);
}

TEST_CASE("Image::mean") {
    CHECK(false);
}

TEST_CASE("Image::limits") {
    CHECK(false);
}

TEST_CASE("Image::get_bounds") {
    CHECK(false);
}

TEST_CASE("Image::set_header") {
    CHECK(false);
}

TEST_CASE("Image::setup_bounds") {
    SECTION("correct_bounds") {
        Matrix data = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
        em::Image image(data);

        em::ObjectBounds2D bounds = image.setup_bounds(1);
        REQUIRE(bounds.size_x() == 6);
        CHECK(bounds[0].min == 1);
        CHECK(bounds[0].max == 4);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 3);
        CHECK(bounds[2].min == 2);
        CHECK(bounds[2].max == 4);
        CHECK(bounds[3].min == 1);
        CHECK(bounds[3].max == 4);
        CHECK(bounds[4].min == 1);
        CHECK(bounds[4].max == 3);
        CHECK(bounds[5].min == 1);
        CHECK(bounds[5].max == 5);

        bounds = image.setup_bounds(5);
        REQUIRE(bounds.size_x() == 6);
        CHECK(bounds[0].min == 3);
        CHECK(bounds[0].max == 3);
        CHECK(bounds[1].min == 2);
        CHECK(bounds[1].max == 3);
        CHECK(bounds[2].min == 0);
        CHECK(bounds[2].max == 0);
        CHECK(bounds[3].min == 3);
        CHECK(bounds[3].max == 3);
        CHECK(bounds[4].min == 3);
        CHECK(bounds[4].max == 3);
        CHECK(bounds[5].min == 5);
        CHECK(bounds[5].max == 5);
    }

    SECTION("more bounds") {
        Matrix data2 = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
        em::Image image2(data2);

        // cutoff = 1
        em::ObjectBounds2D bounds = image2.setup_bounds(1);
        REQUIRE(bounds.size_x() == 6);
        CHECK(bounds[0].min == 1);
        CHECK(bounds[0].max == 5);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 4);
        CHECK(bounds[2].min == 1);
        CHECK(bounds[2].max == 4);
        CHECK(bounds[3].min == 0);
        CHECK(bounds[3].max == 4);
        CHECK(bounds[4].min == 1);
        CHECK(bounds[4].max == 4);
        CHECK(bounds[5].min == 0);
        CHECK(bounds[5].max == 5);

        // cutoff = 2
        bounds = image2.setup_bounds(2);
        REQUIRE(bounds.size_x() == 6);
        CHECK(bounds[0].min == 2);
        CHECK(bounds[0].max == 4);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 4);
        CHECK(bounds[2].min == 2);
        CHECK(bounds[2].max == 2);
        CHECK(bounds[3].min == 0);
        CHECK(bounds[3].max == 3);
        CHECK(bounds[4].min == 2);
        CHECK(bounds[4].max == 2);
        CHECK(bounds[5].min == 0);
        CHECK(bounds[5].max == 3);

        // cutoff = 3
        bounds = image2.setup_bounds(3);
        REQUIRE(bounds.size_x() == 6);
        CHECK(bounds[0].min == 3);
        CHECK(bounds[0].max == 3);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 4);
        CHECK(bounds[2].min == 0);
        CHECK(bounds[2].max == 0);
        CHECK(bounds[3].min == 3);
        CHECK(bounds[3].max == 3);
        CHECK(bounds[4].min == 0);
        CHECK(bounds[4].max == 0);
        CHECK(bounds[5].min == 0);
        CHECK(bounds[5].max == 2);
    }

    SECTION("correct_bounded_area") {
        Matrix data = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
        em::Image image(data);

        em::ObjectBounds2D bounds = image.setup_bounds(1);
        CHECK(bounds.total_area() == 6*6);
        CHECK(bounds.bounded_area() == (4 + 3 + 3 + 4 + 3 + 5));

        bounds = image.setup_bounds(2);
        CHECK(bounds.bounded_area() == (2 + 3 + 2 + 3 + 2 + 3));
    }

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

TEST_CASE("Image::index") {
    CHECK(false);
}

TEST_CASE("Image::squared_sum") {
    CHECK(false);
}

TEST_CASE("Image::operator==") {
    CHECK(false);
}

TEST_CASE("Image::to_string") {
    CHECK(false);
}