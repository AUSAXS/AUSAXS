#include <catch2/catch_test_macros.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <hydrate/generation/RadialHydration.h>
#include <settings/All.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::data;

namespace {
    // a 2x2x2 cube of carbons centred on c
    std::vector<AtomFF> cube(const Vector3<double>& c) {
        std::vector<AtomFF> atoms;
        for (int x : {-1, 1}) {
            for (int y : {-1, 1}) {
                for (int z : {-1, 1}) {
                    atoms.emplace_back(c + Vector3<double>{1.*x, 1.*y, 1.*z}, form_factor::form_factor_t::C);
                }
            }
        }
        return atoms;
    }

    bool all_within(const std::vector<Water>& waters, double distance, const std::vector<Vector3<double>>& centers) {
        return std::all_of(waters.begin(), waters.end(), [&] (const Water& w) {
            return std::any_of(centers.begin(), centers.end(), [&] (const Vector3<double>& c) {
                return w.coordinates().distance(c) < distance;
            });
        });
    }
}

TEST_CASE("GridBasedHydration::hydrate: each body is hydrated from its own atoms") {
    settings::molecule::center = false;
    settings::molecule::implicit_hydrogens = false;
    settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::RadialStrategy;
    hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});

    // body 0 sits at the origin with a single symmetry copy at +60x, body 1 sits at -60x. The grid stores body 0 and its copy as one
    // contiguous block, so a span reconstructed by accumulating size_atom() would hand body 1 the copy's atoms rather than its own.
    Molecule protein({Body{cube({0, 0, 0})}, Body{cube({-60, 0, 0})}});
    protein.get_body(0).symmetry().add(std::make_unique<symmetry::CyclicSymmetry>(
        Vector3<double>{0, 0, 0}, Vector3<double>{60, 0, 0}, Vector3<double>{0, 0, 1}, 0, 1
    ));
    REQUIRE(protein.get_grid()->a_members.size() == 24);

    protein.generate_new_hydration();

    auto w0 = protein.get_body(0).get_waters();
    auto w1 = protein.get_body(1).get_waters();
    REQUIRE(w0.has_value());
    REQUIRE(w1.has_value());
    REQUIRE(0 < w0.value().get().size());
    REQUIRE(0 < w1.value().get().size());

    CHECK(all_within(w0.value().get(), 20, {{0, 0, 0}, {60, 0, 0}}));
    CHECK(all_within(w1.value().get(), 20, {{-60, 0, 0}}));
}
