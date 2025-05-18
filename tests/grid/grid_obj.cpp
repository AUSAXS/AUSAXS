#include <catch2/catch_test_macros.hpp>

#include <utility/Limit3D.h>
#include <grid/Grid.h>

using namespace ausaxs;
using namespace ausaxs::grid;
using namespace ausaxs::grid::detail;

TEST_CASE("GridObj: comparisons") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);
    auto& gref = grid.grid;

    {
        gref.index(0, 0, 0) = A_CENTER;
        CHECK(gref.is_atom_center(0, 0, 0) == true);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == false);

        CHECK(gref.is_only_empty_or_volume(0, 0, 0) == false);
        CHECK(gref.is_only_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_only_atom_center(0, 0, 0) == true);
        CHECK(gref.is_only_volume(0, 0, 0) == false);
    }
    {
        gref.index(0, 0, 0) = A_AREA;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == true);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == true);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == false);

        CHECK(gref.is_only_empty_or_volume(0, 0, 0) == false);
        CHECK(gref.is_only_atom_area_or_volume(0, 0, 0) == true);
        CHECK(gref.is_only_atom_center(0, 0, 0) == false);
        CHECK(gref.is_only_volume(0, 0, 0) == false);
    }
    {
        gref.index(0, 0, 0) = W_CENTER;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_water_center(0, 0, 0) == true);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == false);

        CHECK(gref.is_only_empty_or_volume(0, 0, 0) == false);
        CHECK(gref.is_only_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_only_atom_center(0, 0, 0) == false);
        CHECK(gref.is_only_volume(0, 0, 0) == false);
    }
    {
        gref.index(0, 0, 0) = W_AREA;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == true);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == false);

        CHECK(gref.is_only_empty_or_volume(0, 0, 0) == false);
        CHECK(gref.is_only_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_only_atom_center(0, 0, 0) == false);
        CHECK(gref.is_only_volume(0, 0, 0) == false);
    }
    {
        gref.index(0, 0, 0) = VOLUME;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == true);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == true);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == true);

        CHECK(gref.is_only_empty_or_volume(0, 0, 0) == true);
        CHECK(gref.is_only_atom_area_or_volume(0, 0, 0) == true);
        CHECK(gref.is_only_atom_center(0, 0, 0) == false);
        CHECK(gref.is_only_volume(0, 0, 0) == true);
    }
    {
        gref.index(0, 0, 0) = EMPTY;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == true);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == true);

        CHECK(gref.is_only_empty_or_volume(0, 0, 0) == true);
        CHECK(gref.is_only_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_only_atom_center(0, 0, 0) == false);
        CHECK(gref.is_only_volume(0, 0, 0) == false);
    }
}

TEST_CASE("GridObj: logical operators") {
    State s = A_CENTER;
    CHECK((s & A_CENTER));
    CHECK_FALSE(s & A_AREA);
    CHECK_FALSE(s & W_CENTER);
    CHECK_FALSE(s & W_AREA);
    CHECK_FALSE(s & VOLUME);
    CHECK_FALSE(s & EMPTY);

    s |= A_AREA;
    CHECK(s & A_CENTER);
    CHECK(s & A_AREA);
    CHECK_FALSE(s & W_CENTER);
    CHECK_FALSE(s & W_AREA);
    CHECK_FALSE(s & VOLUME);
    CHECK_FALSE(s & EMPTY);

    s &= ~A_CENTER;
    CHECK_FALSE(s & A_CENTER);
    CHECK(s & A_AREA);
    CHECK_FALSE(s & W_CENTER);
    CHECK_FALSE(s & W_AREA);
    CHECK_FALSE(s & VOLUME);
    CHECK_FALSE(s & EMPTY);
}