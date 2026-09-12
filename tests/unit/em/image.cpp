#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/Image.h>
#include <em/ObjectBounds2D.h>
#include <numeric>
#include <utility/Limit.h>

using namespace ausaxs;

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
    double min = *std::ranges::min_element(data);
    double max = *std::ranges::max_element(data);
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
        CHECK(bounds[0].max == 5);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 4);
        CHECK(bounds[2].min == 2);
        CHECK(bounds[2].max == 5);
        CHECK(bounds[3].min == 1);
        CHECK(bounds[3].max == 5);
        CHECK(bounds[4].min == 1);
        CHECK(bounds[4].max == 4);
        CHECK(bounds[5].min == 1);
        CHECK(bounds[5].max == 6);
        CHECK(image.get_bounds() == bounds);

        bounds = image.setup_bounds(5);
        REQUIRE(bounds.size_x() == 6);
        CHECK(bounds[0].min == 3);
        CHECK(bounds[0].max == 4);
        CHECK(bounds[1].min == 2);
        CHECK(bounds[1].max == 4);
        CHECK(bounds[2].min == 0);
        CHECK(bounds[2].max == 0);   // no voxel above the cutoff, i.e. an empty range
        CHECK(bounds[3].min == 3);
        CHECK(bounds[3].max == 4);
        CHECK(bounds[4].min == 3);
        CHECK(bounds[4].max == 4);
        CHECK(bounds[5].min == 5);
        CHECK(bounds[5].max == 6);
        CHECK(image.get_bounds() == bounds);
    }

    SECTION("more bounds") {
        Matrix data = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
        em::Image image(data);

        em::ObjectBounds2D bounds = image.setup_bounds(1);
        REQUIRE(bounds.size_x() == 6);
        CHECK(bounds[0].min == 1);
        CHECK(bounds[0].max == 6);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 5);
        CHECK(bounds[2].min == 1);
        CHECK(bounds[2].max == 5);
        CHECK(bounds[3].min == 0);
        CHECK(bounds[3].max == 5);
        CHECK(bounds[4].min == 1);
        CHECK(bounds[4].max == 5);
        CHECK(bounds[5].min == 0);
        CHECK(bounds[5].max == 6);
        CHECK(image.get_bounds() == bounds);

        bounds = image.setup_bounds(2);
        REQUIRE(bounds.size_x() == 6);
        CHECK(bounds[0].min == 2);
        CHECK(bounds[0].max == 5);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 5);
        CHECK(bounds[2].min == 2);
        CHECK(bounds[2].max == 3);
        CHECK(bounds[3].min == 0);
        CHECK(bounds[3].max == 4);
        CHECK(bounds[4].min == 2);
        CHECK(bounds[4].max == 3);
        CHECK(bounds[5].min == 0);
        CHECK(bounds[5].max == 4);
        CHECK(image.get_bounds() == bounds);

        bounds = image.setup_bounds(3);
        REQUIRE(bounds.size_x() == 6);
        CHECK(bounds[0].min == 3);
        CHECK(bounds[0].max == 4);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 5);
        CHECK(bounds[2].min == 0);
        CHECK(bounds[2].max == 0);   // no voxel above the cutoff, i.e. an empty range
        CHECK(bounds[3].min == 3);
        CHECK(bounds[3].max == 4);
        CHECK(bounds[4].min == 0);
        CHECK(bounds[4].max == 0);   // no voxel above the cutoff, i.e. an empty range
        CHECK(bounds[5].min == 0);
        CHECK(bounds[5].max == 3);
        CHECK(image.get_bounds() == bounds);
    }

    SECTION("bounds do not change the voxel count") {
        // the bounds are a pure optimisation, so setting them must not change how many voxels are found above the
        // cutoff. At cutoff 5 four of these rows hold a single qualifying voxel and one holds none at all.
        Matrix data = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
        em::Image image(data);

        for (double cutoff : {1., 3., 5.}) {
            int expected = 0;
            for (int x = 0; x < data.N; x++) {
                for (int y = 0; y < data.M; y++) {
                    if (cutoff <= data.index(x, y)) {expected++;}
                }
            }

            image.setup_bounds(cutoff);
            CHECK(image.count_voxels(cutoff) == expected);
        }
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
}

TEST_CASE("Image::index") {
    Matrix<float> data = {
        {1.2, 3.4, 5.6, 7.8, 9.0, 2.3},
        {4.5, 6.7, 8.9, 1.2, 3.4, 5.6},
        {7.8, 9.0, 2.3, 4.5, 6.7, 8.9},
        {1.1, 2.2, 3.3, 4.4, 5.5, 6.6},
        {7.7, 8.8, 9.9, 1.0, 2.2, 3.3},
        {4.4, 5.5, 6.6, 7.7, 8.8, 9.9}
    };
    em::Image image(data);

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
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
