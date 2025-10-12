#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <rigidbody/BodySplitter.h>
#include <grid/Grid.h>
#include <hydrate/generation/RadialHydration.h>
#include <utility/Console.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("RadialHydration: consistency") {
    settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::RadialStrategy;
    hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});

    SECTION("single body") {
        Molecule protein("tests/files/2epe.pdb");
        protein.generate_new_hydration();
        auto h1 = protein.get_waters();
        REQUIRE(h1.size() != 0);
    
        for (int i = 0; i < 10; ++i) {
            protein.generate_new_hydration();
            REQUIRE(h1.size() == protein.size_water());
    
            // check exact equivalence of generated hydration
            auto h2 = protein.get_waters();
            for (unsigned int j = 0; j < h1.size(); ++j) {
                REQUIRE(h1[j].coords == h2[j].coords);
            }
        }    
    }

    SECTION("multiple bodies") {
        Molecule protein = rigidbody::BodySplitter::split("tests/files/2epe.pdb", {20, 40, 60, 80, 100});
        protein.generate_new_hydration();
        auto h1 = protein.get_waters();
        REQUIRE(h1.size() != 0);

        for (int i = 0; i < 10; ++i) {
            protein.generate_new_hydration();
            REQUIRE(h1.size() == protein.size_water());
    
            // check exact equivalence of generated hydration
            auto h2 = protein.get_waters();
            for (unsigned int j = 0; j < h1.size(); ++j) {
                REQUIRE(h1[j].coords == h2[j].coords);
            }
        }
    }
}