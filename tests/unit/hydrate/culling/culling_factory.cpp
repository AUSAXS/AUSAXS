#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <hydrate/culling/BodyCounterCulling.h>
#include <hydrate/culling/CounterCulling.h>
#include <hydrate/culling/CullingFactory.h>
#include <hydrate/culling/NoCulling.h>
#include <hydrate/culling/OutlierCulling.h>
#include <settings/MoleculeSettings.h>

using namespace ausaxs;
using namespace ausaxs::hydrate;
using namespace ausaxs::data;

TEST_CASE("CullingFactory::construct_culling_strategy with global flag") {
    Molecule molecule("tests/files/2epe.pdb");

    SECTION("CounterStrategy global") {
        auto strategy = factory::construct_culling_strategy(
            observer_ptr<Molecule>(&molecule), 
            settings::hydrate::CullingStrategy::CounterStrategy
        );
        REQUIRE(strategy != nullptr);
        auto* counter = dynamic_cast<CounterCulling*>(strategy.get());
        auto* body_counter = dynamic_cast<BodyCounterCulling*>(strategy.get());
        CHECK((counter != nullptr || body_counter != nullptr));
    }

    SECTION("OutlierStrategy") {
        auto strategy = factory::construct_culling_strategy(
            observer_ptr<Molecule>(&molecule), 
            settings::hydrate::CullingStrategy::OutlierStrategy
        );
        REQUIRE(strategy != nullptr);
        auto* outlier = dynamic_cast<OutlierCulling*>(strategy.get());
        CHECK(outlier != nullptr);
    }

    SECTION("NoStrategy") {
        auto strategy = factory::construct_culling_strategy(
            observer_ptr<Molecule>(&molecule), 
            settings::hydrate::CullingStrategy::NoStrategy
        );
        REQUIRE(strategy != nullptr);
        auto* no_culling = dynamic_cast<NoCulling*>(strategy.get());
        CHECK(no_culling != nullptr);
    }
}

TEST_CASE("CullingFactory::construct_culling_strategy with choice") {
    Molecule molecule("tests/files/2epe.pdb");

    SECTION("create with global=true for CounterStrategy") {
        auto strategy = factory::construct_culling_strategy(
            observer_ptr<Molecule>(&molecule), 
            true  // global
        );
        REQUIRE(strategy != nullptr);
    }

    SECTION("create with global=false for CounterStrategy") {
        auto strategy = factory::construct_culling_strategy(
            observer_ptr<Molecule>(&molecule), 
            false  // not global
        );
        REQUIRE(strategy != nullptr);
    }
}
