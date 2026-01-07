#include <catch2/catch_test_macros.hpp>

#include <em/detail/ImageStackBase.h>
#include <em/Image.h>

using namespace ausaxs;

Matrix<float> dummy_image1 = {
    {1, 2, 3},
    {4, 5, 6},
    {7, 8, 9}
};

Matrix<float> dummy_image2 = {
    {10, 11, 12},
    {13, 14, 15},
    {16, 17, 18}
};

Matrix<float> dummy_image3 = {
    {19, 20, 21},
    {22, 23, 24},
    {25, 26, 27}
};

TEST_CASE("ImageStackBase::ImageStackBase") {
    SECTION("std::vector<Image>&") {
        std::vector<em::Image> images;
        for (int i = 0; i < 10; ++i) {
            images.emplace_back(dummy_image1);
        }
        em::ImageStackBase isb(images);
        REQUIRE(isb.size() == 10);
    }
}

TEST_CASE("ImageStackBase::image") {
    std::vector<em::Image> images;
    images.emplace_back(dummy_image1);
    images.emplace_back(dummy_image2);
    images.emplace_back(dummy_image3);
    images[0].set_z(0);
    images[1].set_z(1);
    images[2].set_z(2);

    em::ImageStackBase isb(images);
    REQUIRE(isb.image(0) == images[0]);
    REQUIRE(isb.image(1) == images[1]);
    REQUIRE(isb.image(2) == images[2]);
}

TEST_CASE("ImageStackBase::images") {
    std::vector<em::Image> images;
    images.emplace_back(dummy_image1);
    images.emplace_back(dummy_image2);
    images.emplace_back(dummy_image3);
    images[0].set_z(0);
    images[1].set_z(1);
    images[2].set_z(2);

    em::ImageStackBase isb(images);
    REQUIRE(isb.images() == images);
}

TEST_CASE("ImageStackBase::size") {
    std::vector<em::Image> images;
    images.emplace_back(dummy_image1);
    images.emplace_back(dummy_image2);
    images.emplace_back(dummy_image3);

    em::ImageStackBase isb(images);
    REQUIRE(isb.size() == 3);
}
