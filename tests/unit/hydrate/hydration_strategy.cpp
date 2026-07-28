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
#include <data/symmetry/CyclicSymmetry.h>

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

TEST_CASE("GridBasedHydration::hydrate: each body is hydrated around its own atoms") {
    // The grid stores each body's symmetry replicas alongside its own atoms, so a body's offset into the grid is not the sum of the preceding body sizes.
    // Getting that wrong hydrates a body around some *other* body's replica, which leaves the tail of the molecule bare.
    RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});

    auto cluster = [] (Vector3<double> origin) {
        std::vector<AtomFF> atoms;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    atoms.emplace_back(origin + Vector3<double>{3.0*i, 3.0*j, 3.0*k}, form_factor::form_factor_t::C);
                }
            }
        }
        return atoms;
    };

    // two identically-sized bodies, plus a replica of the first - all three well clear of one another, so a water can be attributed by position alone
    Molecule molecule({Body{cluster({0, 0, 0})}, Body{cluster({60, 0, 0})}});
    molecule.get_body(0).symmetry().add(std::make_unique<symmetry::CyclicSymmetry>(
        Vector3<double>{0, 0, 0}, Vector3<double>{0, 60, 0}, Vector3<double>{0, 0, 1}, 0, 1
    ));
    REQUIRE(molecule.get_body(0).size_atom() == molecule.get_body(1).size_atom());

    RadialHydration hydration{observer_ptr<Molecule>(&molecule), std::make_unique<NoCulling>(observer_ptr<Molecule>(&molecule))};
    hydration.hydrate();

    for (unsigned int i = 0; i < molecule.size_body(); ++i) {
        auto waters = molecule.get_body(i).get_waters();
        REQUIRE(waters.has_value());
        CHECK_FALSE(waters.value().get().empty()); // no body may be left without a shell

        auto cm = molecule.get_body(i).get_cm();    // the body's own cm; the replicas are 60 Å away
        for (const auto& w : waters.value().get()) {
            CHECK(w.coordinates().distance(cm) < 20);
        }
    }
}
