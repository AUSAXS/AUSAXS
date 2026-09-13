// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/detail/header/HeaderFactory.h>
#include <em/detail/header/MRCHeader.h>
#include <em/detail/header/data/MRCData.h>
#include <io/ExistingFile.h>
#include <utility/Axis3D.h>

#include <fstream>
#include <memory>

using namespace ausaxs;

namespace {
    // read only the 1024-byte header, exactly as ImageStackBase does before it reads the voxels
    std::unique_ptr<em::detail::header::IMapHeader> read_header(const io::ExistingFile& file) {
        auto header = em::detail::factory::create_header(file);
        std::ifstream input(file, std::ios::binary);
        REQUIRE(input.is_open());
        input.read(header->get_data_ptr(), header->get_header_size());
        return header;
    }
}

// test.ccp4 and A2M_2020_Q4.ccp4 are two views of the same 485.1 Å cell: the first stores a 3x3x3
// corner of it, the second all 154^3 sampling points. The voxel width follows the sampling rate of
// the cell (mx), not the number of stored columns (nx), so both must decode to 3.15 Å voxels.
TEST_CASE("MRCHeader::get_axes: voxel width is independent of the stored extent") {
    SECTION("sub-volume") {
        auto header = read_header("tests/files/test.ccp4");
        auto axes = header->get_axes();

        CHECK(axes.x.bins == 3);
        CHECK(axes.y.bins == 3);
        CHECK(axes.z.bins == 3);
        CHECK_THAT(axes.x.width(), Catch::Matchers::WithinAbs(3.15, 1e-5));
        CHECK_THAT(axes.y.width(), Catch::Matchers::WithinAbs(3.15, 1e-5));
        CHECK_THAT(axes.z.width(), Catch::Matchers::WithinAbs(3.15, 1e-5));
    }

    SECTION("full cell") {
        auto header = read_header("tests/files/A2M_2020_Q4.ccp4");
        auto axes = header->get_axes();

        CHECK(axes.x.bins == 154);
        CHECK(axes.y.bins == 154);
        CHECK(axes.z.bins == 154);
        CHECK_THAT(axes.x.width(), Catch::Matchers::WithinAbs(3.15, 1e-5));
        CHECK_THAT(axes.y.width(), Catch::Matchers::WithinAbs(3.15, 1e-5));
        CHECK_THAT(axes.z.width(), Catch::Matchers::WithinAbs(3.15, 1e-5));
    }
}

// nx, ny and nz count columns, rows and sections, and (mapc, mapr, maps) names the crystallographic
// axis each of them spans. Attaching a count to the wrong axis is invisible on a cubic map, so this
// uses three different extents and three different voxel widths.
TEST_CASE("MRCHeader::get_axes: stored counts follow the axis order") {
    em::detail::header::MRCData data;
    data.nx = 10;                                               // columns
    data.ny = 20;                                               // rows
    data.nz = 30;                                               // sections
    data.mx = 30;   data.my = 20;   data.mz = 10;               // sampling intervals along x, y, z
    data.cella_x = 60; data.cella_y = 40; data.cella_z = 10;    // so the voxels are 2 x 2 x 1 Å

    SECTION("identity") {
        data.mapc = 1; data.mapr = 2; data.maps = 3;
        auto axes = em::detail::header::MRCHeader(std::move(data)).get_axes();

        CHECK(axes.x.bins == 10);
        CHECK(axes.y.bins == 20);
        CHECK(axes.z.bins == 30);
        CHECK_THAT(axes.x.span(), Catch::Matchers::WithinAbs(20, 1e-5));
        CHECK_THAT(axes.y.span(), Catch::Matchers::WithinAbs(40, 1e-5));
        CHECK_THAT(axes.z.span(), Catch::Matchers::WithinAbs(30, 1e-5));
    }

    // the order carried by every map in tests/files
    SECTION("transposition") {
        data.mapc = 3; data.mapr = 2; data.maps = 1;
        auto axes = em::detail::header::MRCHeader(std::move(data)).get_axes();

        CHECK(axes.x.bins == 30);   // sections
        CHECK(axes.y.bins == 20);   // rows
        CHECK(axes.z.bins == 10);   // columns
        CHECK_THAT(axes.x.span(), Catch::Matchers::WithinAbs(60, 1e-5));
        CHECK_THAT(axes.y.span(), Catch::Matchers::WithinAbs(40, 1e-5));
        CHECK_THAT(axes.z.span(), Catch::Matchers::WithinAbs(10, 1e-5));
    }

    // a cyclic order is not its own inverse, so it separates the mapping from its reverse
    SECTION("cycle") {
        data.mapc = 2; data.mapr = 3; data.maps = 1;
        auto axes = em::detail::header::MRCHeader(std::move(data)).get_axes();

        CHECK(axes.x.bins == 30);   // sections
        CHECK(axes.y.bins == 10);   // columns
        CHECK(axes.z.bins == 20);   // rows
        CHECK_THAT(axes.x.span(), Catch::Matchers::WithinAbs(60, 1e-5));
        CHECK_THAT(axes.y.span(), Catch::Matchers::WithinAbs(20, 1e-5));
        CHECK_THAT(axes.z.span(), Catch::Matchers::WithinAbs(20, 1e-5));
    }
}

// headers constructed in code assign only the cell size and the stored counts, leaving the sampling
// rate and axis order at zero. The other em tests all rely on that decoding as a full, unpermuted cell.
TEST_CASE("MRCHeader::get_axes: unset sampling rate and axis order") {
    em::detail::header::MRCData data;
    data.nx = 6; data.ny = 6; data.nz = 2;
    data.cella_x = 6; data.cella_y = 6; data.cella_z = 2;

    auto axes = em::detail::header::MRCHeader(std::move(data)).get_axes();
    CHECK(axes.x.bins == 6);
    CHECK(axes.y.bins == 6);
    CHECK(axes.z.bins == 2);
    CHECK_THAT(axes.x.width(), Catch::Matchers::WithinAbs(1, 1e-5));
    CHECK_THAT(axes.y.width(), Catch::Matchers::WithinAbs(1, 1e-5));
    CHECK_THAT(axes.z.width(), Catch::Matchers::WithinAbs(1, 1e-5));
}
