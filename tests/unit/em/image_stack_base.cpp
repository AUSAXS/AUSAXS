#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <em/Image.h>
#include <em/detail/ImageStackBase.h>
#include <em/detail/header/MapHeader.h>
#include <em/detail/header/data/MRCData.h>
#include <io/ExistingFile.h>
#include <utility/Axis3D.h>
#include <utility/Exceptions.h>

#include <support/temp_file.h>

#include <array>
#include <filesystem>
#include <fstream>

using namespace ausaxs;

static Matrix<float> dummy_image1 = {
    {1, 2, 3},
    {4, 5, 6},
    {7, 8, 9}
};

static Matrix<float> dummy_image2 = {
    {10, 11, 12},
    {13, 14, 15},
    {16, 17, 18}
};

static Matrix<float> dummy_image3 = {
    {19, 20, 21},
    {22, 23, 24},
    {25, 26, 27}
};

TEST_CASE("ImageStackBase::ImageStackBase") {
    SECTION("std::vector<Image>&") {
        std::vector<em::Image> images;
        images.reserve(10);
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

namespace {
    // a distinct value for every voxel, so a misplaced one cannot coincide with the one belonging there
    float density(int x, int y, int z) {
        return static_cast<float>(100*x + 10*y + z);
    }

    // write the density on an extent[0] x extent[1] x extent[2] grid of 1 Å voxels, stored in the given axis order,
    // preceded by an extended header of ext_size bytes
    void write_map(const io::File& file, const std::array<int, 3>& extent, const std::array<int, 3>& order, int ext_size = 0) {
        em::detail::header::MRCData data;
        data.mode = 2;                                                      // float32
        data.nsymbt = ext_size;
        data.mapc = order[0]; data.mapr = order[1]; data.maps = order[2];
        data.nx = extent[order[0]-1];                                       // the columns span axis order[0], and so on
        data.ny = extent[order[1]-1];
        data.nz = extent[order[2]-1];
        data.mx = extent[0]; data.my = extent[1]; data.mz = extent[2];       // the map spans its full cell
        data.cella_x = static_cast<float>(extent[0]);
        data.cella_y = static_cast<float>(extent[1]);
        data.cella_z = static_cast<float>(extent[2]);

        std::ofstream out(file.path(), std::ios::binary);
        REQUIRE(out.is_open());
        out.write(reinterpret_cast<const char*>(&data), sizeof(data));

        // the extended header is opaque to us, so fill it with something that would be read as a voxel if it were skipped wrongly
        for (int i = 0; i < ext_size; ++i) {
            char byte = static_cast<char>(0xAB);
            out.write(&byte, 1);
        }

        // sections outermost, then rows, then columns, with counter c driving crystallographic axis order[c]
        for (int s = 0; s < data.nz; ++s) {
            for (int r = 0; r < data.ny; ++r) {
                for (int c = 0; c < data.nx; ++c) {
                    std::array<int, 3> coord = {0, 0, 0};
                    coord[order[0]-1] = c;
                    coord[order[1]-1] = r;
                    coord[order[2]-1] = s;

                    float voxel = density(coord[0], coord[1], coord[2]);
                    out.write(reinterpret_cast<const char*>(&voxel), sizeof(voxel));
                }
            }
        }
    }
}

// The axis order is a transparent encoding detail: one density written in any of the six orders must decode to one stack. A cubic map
// cannot show this, since every axis then has the same extent and a count attached to the wrong axis is still the right number.
TEST_CASE("ImageStackBase::read: the axis order is transparent") {
    std::array<int, 3> extent = {2, 3, 4};
    auto order = GENERATE(
        std::array<int, 3>{1, 2, 3}, std::array<int, 3>{2, 1, 3}, std::array<int, 3>{1, 3, 2},   // the identity and the transpositions
        std::array<int, 3>{3, 2, 1}, std::array<int, 3>{2, 3, 1}, std::array<int, 3>{3, 1, 2}    // the last two are the cyclic orders
    );
    INFO("axis order (" << order[0] << ", " << order[1] << ", " << order[2] << ")");

    test::TempFile file(".mrc");
    write_map(file, extent, order);

    em::ImageStackBase isb{io::ExistingFile(file.path())};
    REQUIRE(isb.size() == extent[2]);

    auto axes = isb.get_header()->get_axes();
    REQUIRE(axes.x.bins == extent[0]);
    REQUIRE(axes.y.bins == extent[1]);
    REQUIRE(axes.z.bins == extent[2]);

    for (int x = 0; x < extent[0]; ++x) {
        for (int y = 0; y < extent[1]; ++y) {
            for (int z = 0; z < extent[2]; ++z) {
                INFO("voxel (" << x << ", " << y << ", " << z << ")");
                CHECK(isb.image(z).index(x, y) == density(x, y, z));
            }
        }
    }
}

// an order that is not a permutation of {1, 2, 3} carries no information, so the map is read as if it were stored in the identity
// order. This is the same fallback get_axes() makes, and headers built in code rely on it, since they leave the order at zero.
TEST_CASE("ImageStackBase::read: an unusable axis order falls back to the identity") {
    std::array<int, 3> extent = {2, 3, 4};
    auto stated = GENERATE(std::array<int, 3>{0, 0, 0}, std::array<int, 3>{1, 1, 1}, std::array<int, 3>{1, 2, 4});
    INFO("stated axis order (" << stated[0] << ", " << stated[1] << ", " << stated[2] << ")");

    test::TempFile file(".mrc");
    write_map(file, extent, {1, 2, 3});

    // overwrite the stated order, leaving the voxels laid out in the identity order
    {
        std::fstream out(file.path(), std::ios::binary | std::ios::in | std::ios::out);
        REQUIRE(out.is_open());
        out.seekp(offsetof(em::detail::header::MRCData, mapc));
        out.write(reinterpret_cast<const char*>(stated.data()), sizeof(stated));
    }

    em::ImageStackBase isb{io::ExistingFile(file.path())};
    REQUIRE(isb.size() == extent[2]);

    for (int x = 0; x < extent[0]; ++x) {
        for (int y = 0; y < extent[1]; ++y) {
            for (int z = 0; z < extent[2]; ++z) {
                INFO("voxel (" << x << ", " << y << ", " << z << ")");
                CHECK(isb.image(z).index(x, y) == density(x, y, z));
            }
        }
    }
}

// a crystallographic CCP4 map carries 80 bytes per symmetry operator between the header and the data, and nsymbt states how many.
// That is an encoding detail the density does not know about, so a map must decode identically whether or not it has one.
TEST_CASE("ImageStackBase::read: the extended header is skipped") {
    std::array<int, 3> extent = {2, 3, 4};
    auto ext_size = GENERATE(0, 80, 800);   // no operators, one, and the ten of a typical space group
    INFO("extended header size " << ext_size);

    test::TempFile file(".ccp4");
    write_map(file, extent, {1, 2, 3}, ext_size);

    em::ImageStackBase isb{io::ExistingFile(file.path())};
    REQUIRE(isb.size() == extent[2]);
    REQUIRE(isb.get_header()->get_extended_header_size() == ext_size);

    for (int x = 0; x < extent[0]; ++x) {
        for (int y = 0; y < extent[1]; ++y) {
            for (int z = 0; z < extent[2]; ++z) {
                INFO("voxel (" << x << ", " << y << ", " << z << ")");
                CHECK(isb.image(z).index(x, y) == density(x, y, z));
            }
        }
    }
}

// a file that stops short of what its own header promises must be rejected, not decoded into whatever the unread voxels happen to hold
TEST_CASE("ImageStackBase::read: a truncated file is rejected") {
    std::array<int, 3> extent = {2, 3, 4};

    SECTION("truncated data section") {
        test::TempFile file(".mrc");
        write_map(file, extent, {1, 2, 3}, 80);

        // drop the last four voxels
        std::filesystem::resize_file(file.path(), std::filesystem::file_size(file.path()) - 4*sizeof(float));
        CHECK_THROWS_AS(em::ImageStackBase{io::ExistingFile(file.path())}, except::io_error);
    }

    SECTION("truncated extended header") {
        test::TempFile file(".mrc");
        write_map(file, extent, {1, 2, 3}, 80);

        // keep the header, but cut the extended header in half and lose the data with it
        std::filesystem::resize_file(file.path(), sizeof(em::detail::header::MRCData) + 40);
        CHECK_THROWS_AS(em::ImageStackBase{io::ExistingFile(file.path())}, except::io_error);
    }

    SECTION("truncated header") {
        test::TempFile file(".mrc");
        write_map(file, extent, {1, 2, 3});

        std::filesystem::resize_file(file.path(), sizeof(em::detail::header::MRCData) - 1);
        CHECK_THROWS_AS(em::ImageStackBase{io::ExistingFile(file.path())}, except::io_error);
    }
}
