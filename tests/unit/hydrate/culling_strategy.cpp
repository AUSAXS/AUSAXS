#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hydrate/culling/CullingStrategy.h>
#include <hydrate/culling/NoCulling.h>
#include <hydrate/culling/CounterCulling.h>
#include <hydrate/culling/BodyCounterCulling.h>
#include <hydrate/culling/OutlierCulling.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <grid/detail/GridMember.h>
#include <data/atoms/Water.h>

using namespace ausaxs;
using namespace ausaxs::hydrate;
using namespace ausaxs::data;
using namespace ausaxs::grid;

TEST_CASE("CullingStrategy::set_target_count") {
    Molecule molecule("tests/files/2epe.pdb");
    NoCulling culling{observer_ptr<Molecule>(&molecule)};

    culling.set_target_count(100);
    // The target count is protected, so we can't directly test it
    // but we can verify the method doesn't throw
    CHECK_NOTHROW(culling.set_target_count(200));
}

TEST_CASE("NoCulling::cull") {
    Molecule molecule("tests/files/2epe.pdb");
    NoCulling culling{observer_ptr<Molecule>(&molecule)};

    // Create some water molecules
    std::vector<GridMember<Water>> waters;
    waters.push_back(GridMember<Water>(Water({1, 2, 3}), Vector3<int>{0, 0, 0}));
    waters.push_back(GridMember<Water>(Water({4, 5, 6}), Vector3<int>{1, 1, 1}));
    waters.push_back(GridMember<Water>(Water({7, 8, 9}), Vector3<int>{2, 2, 2}));

    std::span<GridMember<Water>> water_span(waters);
    size_t original_size = water_span.size();

    // NoCulling should not remove any waters
    culling.cull(water_span);
    CHECK(water_span.size() == original_size);
}

TEST_CASE("CounterCulling::construction") {
    Molecule molecule("tests/files/2epe.pdb");
    CounterCulling culling{observer_ptr<Molecule>(&molecule)};

    SECTION("verify construction") {
        // Just verify it constructs successfully
        CHECK(true);
    }
}

TEST_CASE("OutlierCulling::construction") {
    Molecule molecule("tests/files/2epe.pdb");
    OutlierCulling culling{observer_ptr<Molecule>(&molecule)};

    SECTION("verify construction") {
        // Just verify it constructs successfully
        CHECK(true);
    }
}

TEST_CASE("BodyCounterCulling::construction") {
    Molecule molecule("tests/files/2epe.pdb");
    BodyCounterCulling culling{observer_ptr<Molecule>(&molecule)};

    SECTION("verify construction") {
        // Just verify it constructs successfully
        CHECK(true);
    }
}

TEST_CASE("BodyCounterCulling::set_body_ratios") {
    Molecule molecule("tests/files/2epe.pdb");
    BodyCounterCulling culling{observer_ptr<Molecule>(&molecule)};

    std::vector<double> ratios = {0.5, 0.3, 0.2};
    CHECK_NOTHROW(culling.set_body_ratios(ratios));
}
