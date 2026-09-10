// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <grid/Grid.h>
#include <grid/detail/GridObj.h>
#include <hydrate/generation/RadialHydration.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <settings/All.h>

using namespace ausaxs;

namespace {
    /**
     * @brief Recount the bins covered by atoms from scratch.
     *
     * Grid keeps a running tally instead of counting on demand, so comparing the tally against a full sweep is the
     * only way to catch it drifting away from the grid it claims to describe.
     */
    int count_occupied_bins(const grid::Grid& grid) {
        auto bins = grid.get_bins();
        int n = 0;
        for (int i = 0; i < bins.x(); ++i) {
            for (int j = 0; j < bins.y(); ++j) {
                for (int k = 0; k < bins.z(); ++k) {
                    auto state = grid.index(i, j, k);
                    n += (state & (grid::detail::A_CENTER | grid::detail::A_AREA | grid::detail::VOLUME)) != 0;
                }
            }
        }
        return n;
    }

    // The counter only drifts when a body is removed while one of its bins carries both an atom and a water. How
    // often that happens scales with the bin volume, so the wider widths are not decoration: they are what makes the
    // transform loop below reach the faulty state at all.
    void configure(double cell_width) {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::cell_width = cell_width;
        settings::grid::scaling = 0.25;
        settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::RadialStrategy;
        hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});
    }
}

TEST_CASE("Grid: volume counter is consistent across rigidbody moves") {
    double cell_width = GENERATE(1., 2., 3.);
    configure(cell_width);
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

    auto molecule = rigidbody::BodySplitter::split("tests/files/2epe.pdb", {9, 99});
    molecule.generate_new_hydration();
    REQUIRE(molecule.size_water() != 0);

    rigidbody::Rigidbody rigidbody(std::move(molecule));
    REQUIRE(rigidbody.molecule.get_grid()->get_volume_bins() == count_occupied_bins(*rigidbody.molecule.get_grid()));

    // each step moves one body out of the grid and back in, then regenerates the hydration layer, exactly as
    // SimpleController::prepare_step does. Grid::add re-attaches the body's shell unguarded while it still describes
    // the conformation before the move, so waters can land on atoms that have shifted out from under them. At 1 Å
    // the regenerated shell clears those flags again before the next removal, so only the wider widths actually
    // caught the original bug here - the removal invariant below is what covers the default width.
    for (int step = 0; step < 5; ++step) {
        unsigned int ibody = step % rigidbody.molecule.size_body();
        rigidbody.transformer->apply({{1.6, -1.1, 0.8}, {0.09, 0.04, -0.06}}, ibody);
        rigidbody.molecule.generate_new_hydration();

        // apply() may refresh the grid, which reallocates it, so the pointer has to be re-fetched every step
        auto grid = rigidbody.molecule.get_grid();
        INFO("cell_width = " << cell_width << ", step " << step);
        CHECK(grid->get_volume_bins() == count_occupied_bins(*grid));
    }
}

TEST_CASE("Grid: volume counter returns to zero when every body is removed") {
    double cell_width = GENERATE(1., 2., 3.);
    configure(cell_width);

    // loaded rather than split, so the crystallographic waters in the file are kept. Those enter the grid through
    // the same unguarded Grid::add path and are allowed to share a cell with an atom, which reaches the faulty state
    // at every width including the default.
    data::Molecule molecule("tests/files/2epe.pdb");
    REQUIRE(molecule.size_water() != 0);

    auto grid = molecule.get_grid();
    grid->expand_volume();
    int initial = grid->get_volume_bins();
    REQUIRE(initial == count_occupied_bins(*grid));

    for (int cycle = 0; cycle < 3; ++cycle) {
        INFO("cell_width = " << cell_width << ", cycle " << cycle);
        for (auto& body : molecule.get_bodies()) {grid->remove(body);}
        CHECK(grid->get_volume_bins() == 0);

        for (auto& body : molecule.get_bodies()) {grid->add(body);}
        CHECK(grid->get_volume_bins() == count_occupied_bins(*grid));
    }
}
