#include <catch2/catch_test_macros.hpp>

#include <grid/detail/GridObj.h>
#include <math/Vector3.h>

using namespace ausaxs;
using namespace ausaxs::grid::detail;

TEST_CASE("GridObj::constructor") {
    SECTION("default") {
        GridObj grid;
        REQUIRE(grid.size_x() == 0);
        REQUIRE(grid.size_y() == 0);
        REQUIRE(grid.size_z() == 0);
    }

    SECTION("parameterized") {
        GridObj grid(5, 6, 7);
        REQUIRE(grid.size_x() == 5);
        REQUIRE(grid.size_y() == 6);
        REQUIRE(grid.size_z() == 7);
        
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 6; ++j) {
                for (int k = 0; k < 7; ++k) {
                    REQUIRE(grid.index(i, j, k) == EMPTY);
                }
            }
        }
    }
}

TEST_CASE("GridObj::State enum values") {
    SECTION("bit flags are unique") {
        REQUIRE(EMPTY == 0);
        REQUIRE(VOLUME == (1 << 0));
        REQUIRE(VACUUM == (1 << 1));
        REQUIRE(A_CENTER == (1 << 2));
        REQUIRE(A_AREA == (1 << 3));
        REQUIRE(W_CENTER == (1 << 4));
        REQUIRE(W_AREA == (1 << 5));
    }

    SECTION("bitwise OR operations") {
        State combined = VOLUME | A_CENTER;
        REQUIRE((combined & VOLUME) == VOLUME);
        REQUIRE((combined & A_CENTER) == A_CENTER);
        REQUIRE((combined & EMPTY) == EMPTY);
    }

    SECTION("bitwise AND operations") {
        State combined = VOLUME | A_CENTER | A_AREA;
        REQUIRE((combined & (A_CENTER | A_AREA)) != EMPTY);
        REQUIRE((combined & W_CENTER) == EMPTY);
    }
}

TEST_CASE("GridObj::index") {
    GridObj grid(3, 3, 3);
    
    SECTION("Vector3 accessor") {
        grid.index(Vector3<int>(1, 1, 1)) = A_CENTER;
        REQUIRE(grid.index(Vector3<int>(1, 1, 1)) == A_CENTER);
        REQUIRE(grid.index(1, 1, 1) == A_CENTER);
    }

    SECTION("three parameter accessor") {
        grid.index(2, 1, 0) = W_AREA;
        REQUIRE(grid.index(2, 1, 0) == W_AREA);
    }
}

TEST_CASE("GridObj::is_empty") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE(grid.is_empty(EMPTY));
        REQUIRE_FALSE(grid.is_empty(VOLUME));
        REQUIRE_FALSE(grid.is_empty(A_CENTER));
        REQUIRE_FALSE(grid.is_empty(VOLUME | A_CENTER));
    }

    SECTION("coordinate parameters") {
        REQUIRE(grid.is_empty(0, 0, 0));
        grid.index(1, 1, 1) = A_CENTER;
        REQUIRE_FALSE(grid.is_empty(1, 1, 1));
    }
}

TEST_CASE("GridObj::is_volume") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE_FALSE(grid.is_volume(EMPTY));
        REQUIRE(grid.is_volume(VOLUME));
        REQUIRE_FALSE(grid.is_volume(A_CENTER));
        REQUIRE(grid.is_volume(VOLUME | A_CENTER));
    }

    SECTION("coordinate parameters") {
        REQUIRE_FALSE(grid.is_volume(0, 0, 0));
        grid.index(1, 1, 1) = VOLUME;
        REQUIRE(grid.is_volume(1, 1, 1));
        grid.index(2, 2, 2) = VOLUME | A_AREA;
        REQUIRE(grid.is_volume(2, 2, 2));
    }
}

TEST_CASE("GridObj::is_only_volume") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE_FALSE(grid.is_only_volume(EMPTY));
        REQUIRE(grid.is_only_volume(VOLUME));
        REQUIRE_FALSE(grid.is_only_volume(A_CENTER));
        REQUIRE_FALSE(grid.is_only_volume(VOLUME | A_CENTER));
    }

    SECTION("coordinate parameters") {
        grid.index(1, 1, 1) = VOLUME;
        REQUIRE(grid.is_only_volume(1, 1, 1));
        grid.index(2, 2, 2) = VOLUME | A_AREA;
        REQUIRE_FALSE(grid.is_only_volume(2, 2, 2));
    }
}

TEST_CASE("GridObj::is_empty_or_volume") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE(grid.is_empty_or_volume(EMPTY));
        REQUIRE(grid.is_empty_or_volume(VOLUME));
        REQUIRE_FALSE(grid.is_empty_or_volume(A_CENTER));
        REQUIRE(grid.is_empty_or_volume(VOLUME | A_CENTER));
    }

    SECTION("coordinate parameters") {
        REQUIRE(grid.is_empty_or_volume(0, 0, 0));
        grid.index(1, 1, 1) = VOLUME;
        REQUIRE(grid.is_empty_or_volume(1, 1, 1));
        grid.index(2, 2, 2) = A_CENTER;
        REQUIRE_FALSE(grid.is_empty_or_volume(2, 2, 2));
    }
}

TEST_CASE("GridObj::is_only_empty_or_volume") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE(grid.is_only_empty_or_volume(EMPTY));
        REQUIRE(grid.is_only_empty_or_volume(VOLUME));
        REQUIRE_FALSE(grid.is_only_empty_or_volume(A_CENTER));
        REQUIRE_FALSE(grid.is_only_empty_or_volume(VOLUME | A_CENTER));
    }

    SECTION("coordinate parameters") {
        REQUIRE(grid.is_only_empty_or_volume(0, 0, 0));
        grid.index(1, 1, 1) = VOLUME;
        REQUIRE(grid.is_only_empty_or_volume(1, 1, 1));
        grid.index(2, 2, 2) = VOLUME | A_AREA;
        REQUIRE_FALSE(grid.is_only_empty_or_volume(2, 2, 2));
    }
}

TEST_CASE("GridObj::is_empty_or_water") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE(grid.is_empty_or_water(EMPTY));
        REQUIRE_FALSE(grid.is_empty_or_water(VOLUME));
        REQUIRE(grid.is_empty_or_water(W_CENTER));
        REQUIRE(grid.is_empty_or_water(W_AREA));
        REQUIRE(grid.is_empty_or_water(W_CENTER | W_AREA));
        REQUIRE_FALSE(grid.is_empty_or_water(A_CENTER));
    }

    SECTION("coordinate parameters") {
        REQUIRE(grid.is_empty_or_water(0, 0, 0));
        grid.index(1, 1, 1) = W_CENTER;
        REQUIRE(grid.is_empty_or_water(1, 1, 1));
        grid.index(2, 2, 2) = A_CENTER;
        REQUIRE_FALSE(grid.is_empty_or_water(2, 2, 2));
    }
}

TEST_CASE("GridObj::is_empty_or_volume_or_water") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE(grid.is_empty_or_volume_or_water(EMPTY));
        REQUIRE(grid.is_empty_or_volume_or_water(VOLUME));
        REQUIRE(grid.is_empty_or_volume_or_water(W_CENTER));
        REQUIRE(grid.is_empty_or_volume_or_water(W_AREA));
        REQUIRE_FALSE(grid.is_empty_or_volume_or_water(A_CENTER));
        REQUIRE(grid.is_empty_or_volume_or_water(VOLUME | W_CENTER));
    }

    SECTION("coordinate parameters") {
        REQUIRE(grid.is_empty_or_volume_or_water(0, 0, 0));
        grid.index(1, 1, 1) = VOLUME;
        REQUIRE(grid.is_empty_or_volume_or_water(1, 1, 1));
        grid.index(2, 2, 2) = W_AREA;
        REQUIRE(grid.is_empty_or_volume_or_water(2, 2, 2));
    }
}

TEST_CASE("GridObj::is_atom_area") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE_FALSE(grid.is_atom_area(EMPTY));
        REQUIRE(grid.is_atom_area(A_AREA));
        REQUIRE_FALSE(grid.is_atom_area(A_CENTER));
        REQUIRE(grid.is_atom_area(A_AREA | A_CENTER));
    }

    SECTION("coordinate parameters") {
        REQUIRE_FALSE(grid.is_atom_area(0, 0, 0));
        grid.index(1, 1, 1) = A_AREA;
        REQUIRE(grid.is_atom_area(1, 1, 1));
    }
}

TEST_CASE("GridObj::is_atom_area_or_volume") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE_FALSE(grid.is_atom_area_or_volume(EMPTY));
        REQUIRE(grid.is_atom_area_or_volume(A_AREA));
        REQUIRE(grid.is_atom_area_or_volume(VOLUME));
        REQUIRE(grid.is_atom_area_or_volume(A_AREA | VOLUME));
        REQUIRE_FALSE(grid.is_atom_area_or_volume(W_AREA));
    }

    SECTION("coordinate parameters") {
        REQUIRE_FALSE(grid.is_atom_area_or_volume(0, 0, 0));
        grid.index(1, 1, 1) = A_AREA;
        REQUIRE(grid.is_atom_area_or_volume(1, 1, 1));
        grid.index(2, 2, 2) = VOLUME;
        REQUIRE(grid.is_atom_area_or_volume(2, 2, 2));
    }
}

TEST_CASE("GridObj::is_only_atom_area_or_volume") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE_FALSE(grid.is_only_atom_area_or_volume(EMPTY));
        REQUIRE(grid.is_only_atom_area_or_volume(A_AREA));
        REQUIRE(grid.is_only_atom_area_or_volume(VOLUME));
        REQUIRE(grid.is_only_atom_area_or_volume(A_AREA | VOLUME));
        REQUIRE_FALSE(grid.is_only_atom_area_or_volume(A_AREA | W_CENTER));
    }

    SECTION("coordinate parameters") {
        grid.index(1, 1, 1) = A_AREA;
        REQUIRE(grid.is_only_atom_area_or_volume(1, 1, 1));
        grid.index(2, 2, 2) = A_AREA | W_CENTER;
        REQUIRE_FALSE(grid.is_only_atom_area_or_volume(2, 2, 2));
    }
}

TEST_CASE("GridObj::is_water_area") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE_FALSE(grid.is_water_area(EMPTY));
        REQUIRE(grid.is_water_area(W_AREA));
        REQUIRE_FALSE(grid.is_water_area(W_CENTER));
        REQUIRE(grid.is_water_area(W_AREA | W_CENTER));
    }

    SECTION("coordinate parameters") {
        REQUIRE_FALSE(grid.is_water_area(0, 0, 0));
        grid.index(1, 1, 1) = W_AREA;
        REQUIRE(grid.is_water_area(1, 1, 1));
    }
}

TEST_CASE("GridObj::is_atom_center") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE_FALSE(grid.is_atom_center(EMPTY));
        REQUIRE(grid.is_atom_center(A_CENTER));
        REQUIRE_FALSE(grid.is_atom_center(A_AREA));
        REQUIRE(grid.is_atom_center(A_CENTER | A_AREA));
    }

    SECTION("coordinate parameters") {
        REQUIRE_FALSE(grid.is_atom_center(0, 0, 0));
        grid.index(1, 1, 1) = A_CENTER;
        REQUIRE(grid.is_atom_center(1, 1, 1));
    }
}

TEST_CASE("GridObj::is_only_atom_center") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE_FALSE(grid.is_only_atom_center(EMPTY));
        REQUIRE(grid.is_only_atom_center(A_CENTER));
        REQUIRE_FALSE(grid.is_only_atom_center(A_CENTER | A_AREA));
    }

    SECTION("coordinate parameters") {
        grid.index(1, 1, 1) = A_CENTER;
        REQUIRE(grid.is_only_atom_center(1, 1, 1));
        grid.index(2, 2, 2) = A_CENTER | VOLUME;
        REQUIRE_FALSE(grid.is_only_atom_center(2, 2, 2));
    }
}

TEST_CASE("GridObj::is_water_center") {
    GridObj grid(3, 3, 3);

    SECTION("State parameter") {
        REQUIRE_FALSE(grid.is_water_center(EMPTY));
        REQUIRE(grid.is_water_center(W_CENTER));
        REQUIRE_FALSE(grid.is_water_center(W_AREA));
        REQUIRE(grid.is_water_center(W_CENTER | W_AREA));
    }

    SECTION("coordinate parameters") {
        REQUIRE_FALSE(grid.is_water_center(0, 0, 0));
        grid.index(1, 1, 1) = W_CENTER;
        REQUIRE(grid.is_water_center(1, 1, 1));
    }
}
