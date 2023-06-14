#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/detail/ImageStackBase.h>
#include <em/manager/ProteinManager.h>

#include <fstream>

Matrix<float> dummy_data = {
    {1, 2, 3},    {4, 5, 6},    {7, 8, 9},
    {10, 11, 12}, {13, 14, 15}, {16, 17, 18},
    {19, 20, 21}, {22, 23, 24}, {25, 26, 27}
};

struct fixture {
    fixture() {
        for (int i = 0; i < 10; ++i) {
            images.emplace_back(dummy_data);
        }
    }

    std::vector<em::Image> images;
};

TEST_CASE("ImageStackBase::ImageStackBase") {
    SECTION("std::vector<Image>&") {
        std::vector<em::Image> images;
        for (int i = 0; i < 10; ++i) {
            images.emplace_back(dummy_data);
        }
        em::ImageStackBase isb(images);
        REQUIRE(isb.size() == 10);
        REQUIRE(isb.get_header() == nullptr);
        REQUIRE(isb.get_protein_manager() != nullptr);
    }

    SECTION("io::ExistingFile&") {
        io::ExistingFile file("test/files/A2M_2020_Q4.ccp4");
        em::ImageStackBase isb(file);
        REQUIRE(isb.size() == 154);

        auto header = isb.get_header();
        REQUIRE(header->nx == 154);
        REQUIRE(header->ny == 154);
        REQUIRE(header->nz == 154);
        REQUIRE(header->mode == 2);
        REQUIRE(header->mapc == 3);
        REQUIRE(header->mapr == 2);
        REQUIRE(header->maps == 1);

        REQUIRE(isb.get_protein_manager() != nullptr);
    }
}

TEST_CASE_METHOD(fixture, "ImageStackBase::Image") {
        // make one of the images special so we can compare it later
        em::Image special = em::Image({
            {5, 2, 3},    {4, 5, 6},    {7, 8, 9},
            {10, 11, 12}, {13, 14, 15}, {16, 17, 18},
            {19, 20, 21}, {22, 23, 24}, {25, 26, 27}
        });
        images[5] = special;

        em::ImageStackBase isb(images);
        REQUIRE(isb.image(5) == special);
}

TEST_CASE_METHOD(fixture, "ImageStackBase::images") {
    em::ImageStackBase isb(images);
    REQUIRE(isb.images() == images);
}

TEST_CASE_METHOD(fixture, "ImageStackBase::get_histogram") {
    em::ImageStackBase isb(images);
    REQUIRE(isb.get_histogram(5) == isb.get_protein_manager()->get_histogram(5));
}

TEST_CASE("ImageStackBase::count_voxels") {
    SECTION("single image") {
        em::ImageStackBase isb({dummy_data});
        REQUIRE(isb.count_voxels(0) == 27);
        REQUIRE(isb.count_voxels(10) == 17);
        REQUIRE(isb.count_voxels(20) == 7);
        REQUIRE(isb.count_voxels(30) == 0);
    }
}

TEST_CASE_METHOD(fixture, "ImageStackBase::get_protein") {

}
TEST_CASE("ImageStackBase::get_header") {}
TEST_CASE("ImageStackBase::set_header") {}
TEST_CASE("ImageStackBase::size") {}
TEST_CASE("ImageStackBase::save") {}
TEST_CASE("ImageStackBase::get_limits") {}
TEST_CASE("ImageStackBase::mean") {}
TEST_CASE("ImageStackBase::minimum_volume") {}
TEST_CASE("ImageStackBase::from_level") {}
TEST_CASE("ImageStackBase::to_level") {}
TEST_CASE("ImageStackBase::rms") {}
TEST_CASE("ImageStackBase::get_protein_manager") {}
TEST_CASE("ImageStackBase::set_minimum_bounds") {}

TEST_CASE("ImageStackBase::read") {
    // test that the header is read correctly
    SECTION("correct header") {
        std::string file = "test/files/A2M_2020_Q4.ccp4";
        em::ImageStackBase isb(file);

        auto header = isb.get_header();
        REQUIRE(header->nx == 154);
        REQUIRE(header->ny == 154);
        REQUIRE(header->nz == 154);
        REQUIRE(header->mode == 2);
        REQUIRE(header->mapc == 3);
        REQUIRE(header->mapr == 2);
        REQUIRE(header->maps == 1);
    }

    // we want to test that the read function can correctly read map files with different row/column/layer orderings
    SECTION("order") {
        std::string file = "test/files/A2M_2020_Q4.ccp4";
        em::ImageStackBase isb(file);

        auto header = isb.get_header();
        header->nx = 3;
        header->ny = 3;
        header->nz = 3;

        auto save_test_file = [&] (int mapc, int mapr, int maps) {
            header->mapc = mapc;
            header->mapr = mapr;
            header->maps = maps;

            std::vector<Matrix<float>> data = 
            {   {{1, 2, 3},    {4, 5, 6},    {7, 8, 9}},
                {{10, 11, 12}, {13, 14, 15}, {16, 17, 18}},
                {{19, 20, 21}, {22, 23, 24}, {25, 26, 27}}
            };

            std::ofstream output("test/files/test.ccp4", std::ios::binary);
            output.write(reinterpret_cast<char*>(header.get()), sizeof(*header));
            for (auto& m : data) {
                for (auto& v : m) {
                    output.write(reinterpret_cast<char*>(&v), sizeof(v));
                }
            }
        };

        // x = 1, y = 2, z = 3
        save_test_file(1, 2, 3);
        em::ImageStackBase isb1("test/files/test.ccp4");
        REQUIRE(isb1.size() == 3);
        {
            REQUIRE(isb1.image(0).index(0, 0) == 1);
            REQUIRE(isb1.image(0).index(1, 0) == 2);
            REQUIRE(isb1.image(0).index(2, 0) == 3);
            REQUIRE(isb1.image(0).index(0, 1) == 4);
            REQUIRE(isb1.image(0).index(1, 1) == 5);
            REQUIRE(isb1.image(0).index(2, 1) == 6);
            REQUIRE(isb1.image(0).index(0, 2) == 7);
            REQUIRE(isb1.image(0).index(1, 2) == 8);
            REQUIRE(isb1.image(0).index(2, 2) == 9);
        }
        {
            REQUIRE(isb1.image(1).index(0, 0) == 10);
            REQUIRE(isb1.image(1).index(1, 0) == 11);
            REQUIRE(isb1.image(1).index(2, 0) == 12);
            REQUIRE(isb1.image(1).index(0, 1) == 13);
            REQUIRE(isb1.image(1).index(1, 1) == 14);
            REQUIRE(isb1.image(1).index(2, 1) == 15);
            REQUIRE(isb1.image(1).index(0, 2) == 16);
            REQUIRE(isb1.image(1).index(1, 2) == 17);
            REQUIRE(isb1.image(1).index(2, 2) == 18);
        }
        {
            REQUIRE(isb1.image(2).index(0, 0) == 19);
            REQUIRE(isb1.image(2).index(1, 0) == 20);
            REQUIRE(isb1.image(2).index(2, 0) == 21);
            REQUIRE(isb1.image(2).index(0, 1) == 22);
            REQUIRE(isb1.image(2).index(1, 1) == 23);
            REQUIRE(isb1.image(2).index(2, 1) == 24);
            REQUIRE(isb1.image(2).index(0, 2) == 25);
            REQUIRE(isb1.image(2).index(1, 2) == 26);
            REQUIRE(isb1.image(2).index(2, 2) == 27);
        }


        // x = 2, y = 1, z = 3
        save_test_file(2, 1, 3);
        em::ImageStackBase isb2("test/files/test.ccp4");
        REQUIRE(isb2.size() == 3);
        {
            REQUIRE(isb2.image(0).index(0, 0) == 1);
            REQUIRE(isb2.image(0).index(0, 1) == 2);
            REQUIRE(isb2.image(0).index(0, 2) == 3);
            REQUIRE(isb2.image(0).index(1, 0) == 4);
            REQUIRE(isb2.image(0).index(1, 1) == 5);
            REQUIRE(isb2.image(0).index(1, 2) == 6);
            REQUIRE(isb2.image(0).index(2, 0) == 7);
            REQUIRE(isb2.image(0).index(2, 1) == 8);
            REQUIRE(isb2.image(0).index(2, 2) == 9);
        }
        {
            REQUIRE(isb2.image(1).index(0, 0) == 10);
            REQUIRE(isb2.image(1).index(0, 1) == 11);
            REQUIRE(isb2.image(1).index(0, 2) == 12);
            REQUIRE(isb2.image(1).index(1, 0) == 13);
            REQUIRE(isb2.image(1).index(1, 1) == 14);
            REQUIRE(isb2.image(1).index(1, 2) == 15);
            REQUIRE(isb2.image(1).index(2, 0) == 16);
            REQUIRE(isb2.image(1).index(2, 1) == 17);
            REQUIRE(isb2.image(1).index(2, 2) == 18);
        }
        {
            REQUIRE(isb2.image(2).index(0, 0) == 19);
            REQUIRE(isb2.image(2).index(0, 1) == 20);
            REQUIRE(isb2.image(2).index(0, 2) == 21);
            REQUIRE(isb2.image(2).index(1, 0) == 22);
            REQUIRE(isb2.image(2).index(1, 1) == 23);
            REQUIRE(isb2.image(2).index(1, 2) == 24);
            REQUIRE(isb2.image(2).index(2, 0) == 25);
            REQUIRE(isb2.image(2).index(2, 1) == 26);
            REQUIRE(isb2.image(2).index(2, 2) == 27);
        }


        // x = 3, y = 1, z = 2
        save_test_file(3, 1, 2);
        em::ImageStackBase isb3("test/files/test.ccp4");
        REQUIRE(isb3.size() == 3);
        {
            REQUIRE(isb3.image(0).index(0, 0) == 1);
            REQUIRE(isb3.image(0).index(0, 1) == 2);
            REQUIRE(isb3.image(0).index(0, 2) == 3);
            REQUIRE(isb3.image(1).index(0, 0) == 4);
            REQUIRE(isb3.image(1).index(0, 1) == 5);
            REQUIRE(isb3.image(1).index(0, 2) == 6);
            REQUIRE(isb3.image(2).index(0, 0) == 7);
            REQUIRE(isb3.image(2).index(0, 1) == 8);
            REQUIRE(isb3.image(2).index(0, 2) == 9);
        }
        {
            REQUIRE(isb3.image(0).index(1, 0) == 10);
            REQUIRE(isb3.image(0).index(1, 1) == 11);
            REQUIRE(isb3.image(0).index(1, 2) == 12);
            REQUIRE(isb3.image(1).index(1, 0) == 13);
            REQUIRE(isb3.image(1).index(1, 1) == 14);
            REQUIRE(isb3.image(1).index(1, 2) == 15);
            REQUIRE(isb3.image(2).index(1, 0) == 16);
            REQUIRE(isb3.image(2).index(1, 1) == 17);
            REQUIRE(isb3.image(2).index(1, 2) == 18);
        }
        {
            REQUIRE(isb3.image(0).index(2, 0) == 19);
            REQUIRE(isb3.image(0).index(2, 1) == 20);
            REQUIRE(isb3.image(0).index(2, 2) == 21);
            REQUIRE(isb3.image(1).index(2, 0) == 22);
            REQUIRE(isb3.image(1).index(2, 1) == 23);
            REQUIRE(isb3.image(1).index(2, 2) == 24);
            REQUIRE(isb3.image(2).index(2, 0) == 25);
            REQUIRE(isb3.image(2).index(2, 1) == 26);
            REQUIRE(isb3.image(2).index(2, 2) == 27);
        }


        // x = 3, y = 2, z = 1
        save_test_file(3, 2, 1);
        em::ImageStackBase isb4("test/files/test.ccp4");
        REQUIRE(isb4.size() == 3);
        {
            REQUIRE(isb4.image(0).index(0, 0) == 1);
            REQUIRE(isb4.image(1).index(0, 0) == 2);
            REQUIRE(isb4.image(2).index(0, 0) == 3);
            REQUIRE(isb4.image(0).index(0, 1) == 4);
            REQUIRE(isb4.image(1).index(0, 1) == 5);
            REQUIRE(isb4.image(2).index(0, 1) == 6);
            REQUIRE(isb4.image(0).index(0, 2) == 7);
            REQUIRE(isb4.image(1).index(0, 2) == 8);
            REQUIRE(isb4.image(2).index(0, 2) == 9);
        }
        {
            REQUIRE(isb4.image(0).index(1, 0) == 10);
            REQUIRE(isb4.image(1).index(1, 0) == 11);
            REQUIRE(isb4.image(2).index(1, 0) == 12);
            REQUIRE(isb4.image(0).index(1, 1) == 13);
            REQUIRE(isb4.image(1).index(1, 1) == 14);
            REQUIRE(isb4.image(2).index(1, 1) == 15);
            REQUIRE(isb4.image(0).index(1, 2) == 16);
            REQUIRE(isb4.image(1).index(1, 2) == 17);
            REQUIRE(isb4.image(2).index(1, 2) == 18);
        }
        {
            REQUIRE(isb4.image(0).index(2, 0) == 19);
            REQUIRE(isb4.image(1).index(2, 0) == 20);
            REQUIRE(isb4.image(2).index(2, 0) == 21);
            REQUIRE(isb4.image(0).index(2, 1) == 22);
            REQUIRE(isb4.image(1).index(2, 1) == 23);
            REQUIRE(isb4.image(2).index(2, 1) == 24);
            REQUIRE(isb4.image(0).index(2, 2) == 25);
            REQUIRE(isb4.image(1).index(2, 2) == 26);
            REQUIRE(isb4.image(2).index(2, 2) == 27);
        }
    }
}

