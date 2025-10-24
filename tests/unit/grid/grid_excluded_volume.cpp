#include <catch2/catch_test_macros.hpp>

#include <grid/detail/GridExcludedVolume.h>
#include <math/Vector3.h>

using namespace ausaxs;
using namespace ausaxs::grid::exv;

TEST_CASE("GridExcludedVolume::constructor") {
    SECTION("default") {
        GridExcludedVolume exv;
        REQUIRE(exv.interior.empty());
        REQUIRE(exv.surface.empty());
    }
}

TEST_CASE("GridExcludedVolume::has_surface") {
    SECTION("empty") {
        GridExcludedVolume exv;
        REQUIRE_FALSE(exv.has_surface());
    }

    SECTION("with surface") {
        GridExcludedVolume exv;
        exv.surface.push_back({1.0, 2.0, 3.0});
        REQUIRE(exv.has_surface());
    }

    SECTION("with only interior") {
        GridExcludedVolume exv;
        exv.interior.push_back({1.0, 2.0, 3.0});
        REQUIRE_FALSE(exv.has_surface());
    }
}

TEST_CASE("GridExcludedVolume::vectors") {
    GridExcludedVolume exv;

    SECTION("interior") {
        exv.interior.push_back({1.0, 2.0, 3.0});
        exv.interior.push_back({4.0, 5.0, 6.0});
        
        REQUIRE(exv.interior.size() == 2);
        REQUIRE(exv.interior[0] == Vector3<double>(1.0, 2.0, 3.0));
        REQUIRE(exv.interior[1] == Vector3<double>(4.0, 5.0, 6.0));
    }

    SECTION("surface") {
        exv.surface.push_back({1.0, 2.0, 3.0});
        exv.surface.push_back({4.0, 5.0, 6.0});
        
        REQUIRE(exv.surface.size() == 2);
        REQUIRE(exv.surface[0] == Vector3<double>(1.0, 2.0, 3.0));
        REQUIRE(exv.surface[1] == Vector3<double>(4.0, 5.0, 6.0));
    }
}
