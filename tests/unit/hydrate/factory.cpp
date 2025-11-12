#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hydrate/culling/CullingFactory.h>
#include <hydrate/culling/CounterCulling.h>
#include <hydrate/culling/BodyCounterCulling.h>
#include <hydrate/culling/OutlierCulling.h>
#include <hydrate/culling/NoCulling.h>
#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/generation/RadialHydration.h>
#include <hydrate/generation/AxesHydration.h>
#include <hydrate/generation/PepsiHydration.h>
#include <hydrate/generation/JanHydration.h>
#include <hydrate/generation/NoHydration.h>
#include <data/Molecule.h>
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

TEST_CASE("HydrationFactory::construct_hydration_generator") {
    Molecule molecule("tests/files/2epe.pdb");

    SECTION("AxesStrategy") {
        auto generator = factory::construct_hydration_generator(
            observer_ptr<Molecule>(&molecule),
            settings::hydrate::HydrationStrategy::AxesStrategy
        );
        REQUIRE(generator != nullptr);
        auto* axes = dynamic_cast<AxesHydration*>(generator.get());
        CHECK(axes != nullptr);
    }

    SECTION("RadialStrategy") {
        auto generator = factory::construct_hydration_generator(
            observer_ptr<Molecule>(&molecule),
            settings::hydrate::HydrationStrategy::RadialStrategy
        );
        REQUIRE(generator != nullptr);
        auto* radial = dynamic_cast<RadialHydration*>(generator.get());
        CHECK(radial != nullptr);
    }

    SECTION("JanStrategy") {
        auto generator = factory::construct_hydration_generator(
            observer_ptr<Molecule>(&molecule),
            settings::hydrate::HydrationStrategy::JanStrategy
        );
        REQUIRE(generator != nullptr);
        auto* jan = dynamic_cast<JanHydration*>(generator.get());
        CHECK(jan != nullptr);
    }

    SECTION("NoStrategy") {
        auto generator = factory::construct_hydration_generator(
            observer_ptr<Molecule>(&molecule),
            settings::hydrate::HydrationStrategy::NoStrategy
        );
        REQUIRE(generator != nullptr);
        auto* no_hydration = dynamic_cast<NoHydration*>(generator.get());
        CHECK(no_hydration != nullptr);
    }

    SECTION("PepsiStrategy") {
        auto generator = factory::construct_hydration_generator(
            observer_ptr<Molecule>(&molecule),
            settings::hydrate::HydrationStrategy::PepsiStrategy
        );
        REQUIRE(generator != nullptr);
        auto* pepsi = dynamic_cast<PepsiHydration*>(generator.get());
        CHECK(pepsi != nullptr);
    }
}

TEST_CASE("HydrationFactory::construct_hydration_generator with culling strategy") {
    Molecule molecule("tests/files/2epe.pdb");

    SECTION("RadialStrategy with CounterStrategy") {
        auto generator = factory::construct_hydration_generator(
            observer_ptr<Molecule>(&molecule),
            settings::hydrate::HydrationStrategy::RadialStrategy,
            settings::hydrate::CullingStrategy::CounterStrategy
        );
        REQUIRE(generator != nullptr);
        auto* radial = dynamic_cast<RadialHydration*>(generator.get());
        CHECK(radial != nullptr);
    }

    SECTION("AxesStrategy with OutlierStrategy") {
        auto generator = factory::construct_hydration_generator(
            observer_ptr<Molecule>(&molecule),
            settings::hydrate::HydrationStrategy::AxesStrategy,
            settings::hydrate::CullingStrategy::OutlierStrategy
        );
        REQUIRE(generator != nullptr);
        auto* axes = dynamic_cast<AxesHydration*>(generator.get());
        CHECK(axes != nullptr);
    }
}

TEST_CASE("HydrationFactory::construct_hydration_generator with default settings") {
    Molecule molecule("tests/files/2epe.pdb");

    SECTION("use default settings") {
        auto generator = factory::construct_hydration_generator(
            observer_ptr<Molecule>(&molecule)
        );
        REQUIRE(generator != nullptr);
    }
}
