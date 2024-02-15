#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <em/detail/ImageStackBase.h>
#include <em/detail/header/data/MRCData.h>
#include <em/detail/header/MRCHeader.h>
#include <em/manager/ProteinManager.h>
#include <hist/HistFwd.h>
#include <em/ObjectBounds3D.h>
#include <data/Molecule.h>
#include <settings/All.h>

#include <fstream>

Matrix<float> dummy_image_stack = {
    {1, 2, 3},    {4, 5, 6},    {7, 8, 9},
    {10, 11, 12}, {13, 14, 15}, {16, 17, 18},
    {19, 20, 21}, {22, 23, 24}, {25, 26, 27}
};

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

struct fixture {
    fixture() {
        images.emplace_back(dummy_image1);
        images.emplace_back(dummy_image2);
        images.emplace_back(dummy_image3);
        images[0].set_z(0);
        images[1].set_z(1);
        images[2].set_z(2);
    }

    std::vector<em::Image> images;
};

TEST_CASE("ImageStackBase::ImageStackBase") {
    SECTION("std::vector<Image>&") {
        std::vector<em::Image> images;
        for (int i = 0; i < 10; ++i) {
            images.emplace_back(dummy_image1);
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

        auto header = static_cast<em::detail::header::MRCData*>(isb.get_header()->get_data());
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

TEST_CASE_METHOD(fixture, "ImageStackBase::image") {
    em::ImageStackBase isb(images);
    REQUIRE(isb.image(0) == images[0]);
    REQUIRE(isb.image(1) == images[1]);
    REQUIRE(isb.image(2) == images[2]);
}

TEST_CASE_METHOD(fixture, "ImageStackBase::images") {
    em::ImageStackBase isb(images);
    REQUIRE(isb.images() == images);
}

TEST_CASE_METHOD(fixture, "ImageStackBase::get_histogram") {
    settings::molecule::use_effective_charge = false;
    em::ImageStackBase isb("test/files/A2M_2020_Q4.ccp4");
    REQUIRE(isb.get_histogram(5)->get_total_counts() == isb.get_protein_manager()->get_histogram(5)->get_total_counts());
}

TEST_CASE_METHOD(fixture, "ImageStackBase::count_voxels") {
    SECTION("single image") {
        em::ImageStackBase isb(images);
        REQUIRE(isb.count_voxels(0) == 27);
        REQUIRE(isb.count_voxels(10) == 18);
        REQUIRE(isb.count_voxels(20) == 8);
        REQUIRE(isb.count_voxels(30) == 0);
    }
}

TEST_CASE_METHOD(fixture, "ImageStackBase::get_protein") {
    em::ImageStackBase isb("test/files/A2M_2020_Q4.ccp4");
    REQUIRE(isb.get_protein(5) == isb.get_protein_manager()->get_protein(5));
}

TEST_CASE_METHOD(fixture, "ImageStackBase::get_header") {
    SECTION("no header") {
        em::ImageStackBase isb(images);
        REQUIRE(isb.get_header() == nullptr);
    }

    SECTION("with header") {
        io::ExistingFile file("test/files/A2M_2020_Q4.ccp4");
        em::ImageStackBase isb(file);
        REQUIRE(isb.get_header() != nullptr);
    }
}

TEST_CASE_METHOD(fixture, "ImageStackBase::set_header") {
    std::unique_ptr<em::detail::header::MapHeader> header; 
    {
        em::detail::header::MRCData header_data;
        header_data.nx = 10;
        header_data.ny = 10;
        header_data.nz = 10;
        header_data.mode = 2;
        header_data.mapc = 3;
        header_data.mapr = 2;
        header_data.maps = 1;
        header = std::make_unique<em::detail::header::MRCHeader>(std::move(header_data));
    }

    em::ImageStackBase isb(images);
    isb.set_header(std::move(header));
    auto h = static_cast<em::detail::header::MRCData*>(isb.get_header()->get_data());
    CHECK(h->nx == 10);
    CHECK(h->ny == 10);
    CHECK(h->nz == 10);
    CHECK(h->mode == 2);
    CHECK(h->mapc == 3);
    CHECK(h->mapr == 2);
    CHECK(h->maps == 1);
}

TEST_CASE_METHOD(fixture, "ImageStackBase::size") {
    SECTION("single image") {
        em::ImageStackBase isb({dummy_image1});
        REQUIRE(isb.size() == 1);
    }

    SECTION("multiple images") {
        em::ImageStackBase isb(images);
        REQUIRE(isb.size() == 3);
    }
}

TEST_CASE_METHOD(fixture, "ImageStackBase::save") {
    io::File file("test/temp/ImageStackBase.save.pdb");
    em::ImageStackBase isb(images);

    em::detail::header::MRCData header_data;
    header_data.cella_x = 1; header_data.cella_y = 1; header_data.cella_z = 1;
    header_data.nx = 3; header_data.ny = 3; header_data.nz = 3;
    std::unique_ptr header = std::make_unique<em::detail::header::MRCHeader>(std::move(header_data));

    isb.set_header(std::move(header));
    isb.save(5, file);
    REQUIRE(file.exists());
}

TEST_CASE("ImageStackBase::mean") {
    SECTION("single image") {
        em::ImageStackBase isb({dummy_image1});
        REQUIRE(isb.mean() == 5);
    }

    SECTION("multiple images") {
        em::ImageStackBase isb({dummy_image1, dummy_image2, dummy_image3});
        REQUIRE(isb.mean() == 14);
    }
}

TEST_CASE_METHOD(fixture, "ImageStackBase::minimum_volume") {
    em::ImageStackBase isb(images);

    auto val = GENERATE(1, 5, 9);
    auto vol = isb.minimum_volume(val);
    REQUIRE(vol.size_x() == 3);
    REQUIRE(vol[0] == isb.image(0).setup_bounds(val));
    REQUIRE(vol[1] == isb.image(1).setup_bounds(val));
    REQUIRE(vol[2] == isb.image(2).setup_bounds(val));
}

TEST_CASE_METHOD(fixture, "ImageStackBase::from_level") {
    em::ImageStackBase isb(images);
    auto sigma = GENERATE(0.5, 1, 1.5, 2);
    REQUIRE(isb.from_level(sigma) == isb.rms()*sigma);
}

TEST_CASE_METHOD(fixture, "ImageStackBase::to_level") {
    em::ImageStackBase isb(images);
    auto sigma = GENERATE(0.5, 1, 1.5, 2);
    REQUIRE(isb.to_level(sigma) == sigma/isb.rms());
}

TEST_CASE("ImageStackBase::rms") {
    SECTION("single image") {
        em::ImageStackBase isb({dummy_image1});
        REQUIRE_THAT(isb.rms(), Catch::Matchers::WithinAbs(5.62731434, 1e-3));
    }

    SECTION("multiple images") {
        em::ImageStackBase isb({dummy_image1, dummy_image2, dummy_image3});
        REQUIRE_THAT(isb.rms(), Catch::Matchers::WithinAbs(16.0208198, 1e-3));
    }
}

TEST_CASE_METHOD(fixture, "ImageStackBase::set_minimum_bounds") {
    em::ImageStackBase isb(images);

    auto bound = GENERATE(1, 5, 9);
    isb.set_minimum_bounds(bound);
    for (unsigned int i = 0; i < isb.size(); ++i) {
        auto b = isb.image(i).get_bounds();
        for (auto j = b[i].min; j < b[i].max; ++j) {
            REQUIRE(bound <= isb.image(i).index(i, j));
        }
    }
}

TEST_CASE("ImageStackBase::read") {
    // test that the header is read correctly
    SECTION("correct header") {
        std::string file = "test/files/A2M_2020_Q4.ccp4";
        em::ImageStackBase isb(file);

        auto header = static_cast<em::detail::header::MRCData*>(isb.get_header()->get_data());
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

        auto header = static_cast<em::detail::header::MRCData*>(isb.get_header()->get_data());
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
            output.write(reinterpret_cast<char*>(header), sizeof(*header));
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

