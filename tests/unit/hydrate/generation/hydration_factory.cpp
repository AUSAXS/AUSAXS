#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <hydrate/generation/AxesHydration.h>
#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/generation/JanHydration.h>
#include <hydrate/generation/NoHydration.h>
#include <hydrate/generation/PepsiHydration.h>
#include <hydrate/generation/RadialHydration.h>
#include <settings/MoleculeSettings.h>

using namespace ausaxs;
using namespace ausaxs::hydrate;
using namespace ausaxs::data;

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
