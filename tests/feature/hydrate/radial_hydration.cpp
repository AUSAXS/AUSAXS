#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/atoms/Atom.h>
#include <hydrate/generation/RadialHydration.h>
#include <rigidbody/BodySplitter.h>
#include <settings/All.h>

#include <algorithm>
#include <cmath>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("RadialHydration: consistency") {
    settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::RadialStrategy;
    hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});

    SECTION("single body") {
        Molecule protein("tests/files/2epe.pdb");
        protein.generate_new_hydration();
        auto h1 = protein.get_waters();
        REQUIRE(!h1.empty());
    
        for (int i = 0; i < 10; ++i) {
            protein.generate_new_hydration();
            REQUIRE(static_cast<int>(h1.size()) == protein.size_water());
    
            // check exact equivalence of generated hydration
            auto h2 = protein.get_waters();
            for (int j = 0; j < static_cast<int>(h1.size()); ++j) {
                REQUIRE(h1[j].coords == h2[j].coords);
            }
        }    
    }

    SECTION("multiple bodies") {
        Molecule protein = rigidbody::BodySplitter::split("tests/files/2epe.pdb", {20, 40, 60, 80, 100});
        protein.generate_new_hydration();
        auto h1 = protein.get_waters();
        REQUIRE(!h1.empty());

        for (int i = 0; i < 10; ++i) {
            protein.generate_new_hydration();
            REQUIRE(static_cast<int>(h1.size()) == protein.size_water());
    
            // check exact equivalence of generated hydration
            auto h2 = protein.get_waters();
            for (int j = 0; j < static_cast<int>(h1.size()); ++j) {
                REQUIRE(h1[j].coords == h2[j].coords);
            }
        }
    }
}
TEST_CASE("RadialHydration: atoms without a form factor still exclude water") {
    // molecule_from_arrays and the C++ Body(std::vector<Atom>) constructor give every atom the UNKNOWN form factor;
    // the grid must still give those atoms a volume, or the shell is placed inside the molecule
    settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::RadialStrategy;
    hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});

    Molecule file("tests/files/2epe.pdb");
    file.clear_hydration();
    std::vector<Atom> bare;
    for (const auto& a : file.iterate_atoms()) {bare.emplace_back(a.coordinates(), a.weight());}
    Molecule arrays({Body{bare}});

    file.generate_new_hydration();
    arrays.generate_new_hydration();
    REQUIRE(arrays.size_water() != 0);
    CHECK(std::abs(arrays.size_water() - file.size_water()) < 0.25*file.size_water());

    double closest = 1e9;
    for (const auto& w : arrays.get_waters()) {
        for (const auto& a : arrays.iterate_atoms()) {closest = std::min(closest, a.coordinates().distance(w.coordinates()));}
    }
    CHECK(closest > 2);
}
