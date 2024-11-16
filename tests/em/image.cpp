#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/Image.h>
#include <em/ImageStack.h>
#include <em/ObjectBounds3D.h>
#include <numeric>
#include <utility/Limit.h>

using namespace ausaxs;

// TODO: fix the remaining tests
TEST_CASE("Image::Image") {
    SECTION("std::shared_ptr<ccp4::Header>, unsigned int") {}
    SECTION("Matrix<float>&") {}
    SECTION("Matrix<float>&, std::shared_ptr<ccp4::Header>, unsigned int") {}
}

TEST_CASE("Image::as_hist") {}
TEST_CASE("Image::generate_atoms") {}
TEST_CASE("Image::set_header") {}

TEST_CASE("Image::count_voxels") {
    Matrix<float> data = {
        {1, 1, 1, 1, 1, 1}, 
        {2, 2, 2, 2, 2, 2}, 
        {3, 3, 3, 3, 3, 3}, 
        {4, 4, 4, 4, 4, 4}, 
        {5, 5, 5, 5, 5, 5}, 
        {6, 6, 6, 6, 6, 6}
    };

    em::Image image(data);
    CHECK(image.count_voxels(6) == 6);
    CHECK(image.count_voxels(5) == 12);
    CHECK(image.count_voxels(4) == 18);
    CHECK(image.count_voxels(3) == 24);
    CHECK(image.count_voxels(2) == 30);
    CHECK(image.count_voxels(1) == 36);
}

TEST_CASE("Image: get & set_z") {
    em::Image image(Matrix<float>(0, 0));
    image.set_z(5);
    CHECK(image.get_z() == 5);

    image.set_z(10);
    CHECK(image.get_z() == 10);
}

TEST_CASE("Image::mean") {
    Matrix<float> data = {
        {1.2, 3.4, 5.6, 7.8, 9.0, 2.3},
        {4.5, 6.7, 8.9, 1.2, 3.4, 5.6},
        {7.8, 9.0, 2.3, 4.5, 6.7, 8.9},
        {1.1, 2.2, 3.3, 4.4, 5.5, 6.6},
        {7.7, 8.8, 9.9, 1.0, 2.2, 3.3},
        {4.4, 5.5, 6.6, 7.7, 8.8, 9.9}
    };

    em::Image image(data);
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    CHECK_THAT(image.mean(), Catch::Matchers::WithinAbs(sum/(6*6), 1e-3));
}

TEST_CASE("Image::limits") {
    Matrix<float> data = {
        {1.2, 3.4, 5.6, 7.8, 9.0, 2.3},
        {4.5, 6.7, 8.9, 1.2, 3.4, 5.6},
        {7.8, 9.0, 2.3, 4.5, 6.7, 8.9},
        {1.1, 2.2, 3.3, 4.4, 5.5, 6.6},
        {7.7, 8.8, 9.9, 1.0, 2.2, 3.3},
        {4.4, 5.5, 6.6, 7.7, 8.8, 9.9}
    };

    em::Image image(data);
    double min = *std::min_element(data.begin(), data.end());
    double max = *std::max_element(data.begin(), data.end());
    CHECK(image.limits().min == min);
    CHECK(image.limits().max == max);
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
        CHECK(image.get_bounds() == bounds);

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
        CHECK(image.get_bounds() == bounds);
    }

    SECTION("more bounds") {
        Matrix data = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
        em::Image image(data);

        // cutoff = 1
        em::ObjectBounds2D bounds = image.setup_bounds(1);
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
        CHECK(image.get_bounds() == bounds);

        // cutoff = 2
        bounds = image.setup_bounds(2);
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
        CHECK(image.get_bounds() == bounds);

        // cutoff = 3
        bounds = image.setup_bounds(3);
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
        CHECK(image.get_bounds() == bounds);
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
    // use some more creative data 
    Matrix<float> data = {
        {1.2, 3.4, 5.6, 7.8, 9.0, 2.3},
        {4.5, 6.7, 8.9, 1.2, 3.4, 5.6},
        {7.8, 9.0, 2.3, 4.5, 6.7, 8.9},
        {1.1, 2.2, 3.3, 4.4, 5.5, 6.6},
        {7.7, 8.8, 9.9, 1.0, 2.2, 3.3},
        {4.4, 5.5, 6.6, 7.7, 8.8, 9.9}
    };
    em::Image image(data);

    for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            CHECK(image.index(i, j) == data.index(i, j));
        }
    }
}

TEST_CASE("Image::squared_sum") {
    Matrix<float> data = {
        {1.2, 3.4, 5.6, 7.8, 9.0, 2.3},
        {4.5, 6.7, 8.9, 1.2, 3.4, 5.6},
        {7.8, 9.0, 2.3, 4.5, 6.7, 8.9},
        {1.1, 2.2, 3.3, 4.4, 5.5, 6.6},
        {7.7, 8.8, 9.9, 1.0, 2.2, 3.3},
        {4.4, 5.5, 6.6, 7.7, 8.8, 9.9}
    };

    em::Image image(data);
    double sqsum = std::accumulate(data.begin(), data.end(), 0.0, [](double sum, float val) {return sum + val*val;});
    CHECK_THAT(image.squared_sum(), Catch::Matchers::WithinAbs(sqsum, 1e-3));
}

TEST_CASE("Image::operator==") {
    Matrix<float> data1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Matrix<float> data2 = {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};

    em::Image image1(data1);
    em::Image image2(data2);

    CHECK(image1 != image2);

    image2 = image1;
    CHECK(image1 == image2);
}