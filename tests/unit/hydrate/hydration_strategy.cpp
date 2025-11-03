#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hydrate/generation/HydrationStrategy.h>
#include <hydrate/generation/NoHydration.h>
#include <hydrate/generation/GridBasedHydration.h>
#include <hydrate/generation/RadialHydration.h>
#include <hydrate/generation/AxesHydration.h>
#include <hydrate/generation/PepsiHydration.h>
#include <hydrate/generation/JanHydration.h>
#include <hydrate/culling/NoCulling.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::hydrate;
using namespace ausaxs::data;

TEST_CASE("NoHydration::global") {
    NoHydration hydration;
    CHECK(hydration.global() == false);
}

TEST_CASE("NoHydration::hydrate") {
    NoHydration hydration;
    // Should not throw and does nothing
    CHECK_NOTHROW(hydration.hydrate());
}

TEST_CASE("RadialHydration::constructor") {
    Molecule molecule("tests/files/2epe.pdb");

    SECTION("constructor with molecule") {
        RadialHydration hydration{observer_ptr<Molecule>(&molecule)};
        CHECK(hydration.global() == false);
    }

    SECTION("constructor with molecule and culling strategy") {
        auto culling = std::make_unique<NoCulling>(observer_ptr<Molecule>(&molecule));
        RadialHydration hydration{observer_ptr<Molecule>(&molecule), std::move(culling)};
        CHECK(hydration.global() == false);
    }
}

TEST_CASE("RadialHydration::global") {
    Molecule molecule("tests/files/2epe.pdb");
    RadialHydration hydration{observer_ptr<Molecule>(&molecule)};
    CHECK(hydration.global() == false);
}

TEST_CASE("RadialHydration::set_noise_generator") {
    // Test that we can set a custom noise generator
    CHECK_NOTHROW(RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};}));
    CHECK_NOTHROW(RadialHydration::set_noise_generator([] () {return Vector3<double>{1, 1, 1};}));
}

TEST_CASE("AxesHydration::constructor") {
    Molecule molecule("tests/files/2epe.pdb");

    SECTION("constructor with molecule") {
        AxesHydration hydration{observer_ptr<Molecule>(&molecule)};
        CHECK(hydration.global() == false);
    }

    SECTION("constructor with molecule and culling strategy") {
        auto culling = std::make_unique<NoCulling>(observer_ptr<Molecule>(&molecule));
        AxesHydration hydration{observer_ptr<Molecule>(&molecule), std::move(culling)};
        CHECK(hydration.global() == false);
    }
}

TEST_CASE("AxesHydration::global") {
    Molecule molecule("tests/files/2epe.pdb");
    AxesHydration hydration{observer_ptr<Molecule>(&molecule)};
    CHECK(hydration.global() == false);
}

TEST_CASE("PepsiHydration::constructor") {
    Molecule molecule("tests/files/2epe.pdb");

    SECTION("constructor with molecule") {
        PepsiHydration hydration{observer_ptr<Molecule>(&molecule)};
        CHECK(hydration.global() == false);
    }

    SECTION("constructor with molecule and culling strategy") {
        auto culling = std::make_unique<NoCulling>(observer_ptr<Molecule>(&molecule));
        PepsiHydration hydration{observer_ptr<Molecule>(&molecule), std::move(culling)};
        CHECK(hydration.global() == false);
    }
}

TEST_CASE("PepsiHydration::global") {
    Molecule molecule("tests/files/2epe.pdb");
    PepsiHydration hydration{observer_ptr<Molecule>(&molecule)};
    CHECK(hydration.global() == false);
}

TEST_CASE("JanHydration::constructor") {
    Molecule molecule("tests/files/2epe.pdb");

    SECTION("constructor with molecule") {
        JanHydration hydration{observer_ptr<Molecule>(&molecule)};
        CHECK(hydration.global() == true);
    }

    SECTION("constructor with molecule and culling strategy") {
        auto culling = std::make_unique<NoCulling>(observer_ptr<Molecule>(&molecule));
        JanHydration hydration{observer_ptr<Molecule>(&molecule), std::move(culling)};
        CHECK(hydration.global() == true);
    }
}

TEST_CASE("JanHydration::global") {
    Molecule molecule("tests/files/2epe.pdb");
    JanHydration hydration{observer_ptr<Molecule>(&molecule)};
    CHECK(hydration.global() == true);
}

TEST_CASE("GridBasedHydration::set_culling_strategy") {
    Molecule molecule("tests/files/2epe.pdb");
    RadialHydration hydration{observer_ptr<Molecule>(&molecule)};

    auto culling = std::make_unique<NoCulling>(observer_ptr<Molecule>(&molecule));
    CHECK_NOTHROW(hydration.set_culling_strategy(std::move(culling)));
}
