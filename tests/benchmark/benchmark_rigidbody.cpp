// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/controller/SimpleController.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/selection/BodySelectStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/ExistingFile.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;

// Helper: build a fresh Rigidbody + controller from the LAR1-2 test structure (2 bodies).
static Rigidbody make_lar12() {
    Rigidbody rb = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    rb.molecule.generate_new_hydration();
    return rb;
}

// The fixture runs once per BENCHMARK_ADVANCED invocation via BENCHMARK_ADVANCED's setup block.
TEST_CASE("Rigidbody pipeline step benchmarks: LAR1-2 (2 bodies)", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::general::warnings = false;
    settings::molecule::allow_unknown_residues = true;
    settings::molecule::implicit_hydrogens = false;

    // -----------------------------------------------------------------------
    // Full iteration (baseline) — prepare_step + finish_step together.
    // -----------------------------------------------------------------------
    BENCHMARK_ADVANCED("Full iteration (prepare + finish)") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = make_lar12();
        rb.controller->setup(io::ExistingFile("tests/files/LAR1-2.dat"));
        meter.measure([&] {
            rb.controller->prepare_step();
            rb.controller->finish_step();
        });
    };

    // -----------------------------------------------------------------------
    // Body selection.
    // -----------------------------------------------------------------------
    BENCHMARK_ADVANCED("Body selection") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = make_lar12();
        meter.measure([&] {
            return rb.body_selector->next_mask();
        });
    };

    // -----------------------------------------------------------------------
    // Parameter generation.
    // -----------------------------------------------------------------------
    BENCHMARK_ADVANCED("Parameter generation") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = make_lar12();
        auto [ibody, iconstraint, mask] = rb.body_selector->next_mask();
        meter.measure([&] {
            return rb.parameter_generator->next(ibody);
        });
    };

    // -----------------------------------------------------------------------
    // Transform application (body move).
    // -----------------------------------------------------------------------
    BENCHMARK_ADVANCED("Transform apply") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = make_lar12();
        auto [ibody, iconstraint, mask] = rb.body_selector->next_mask();
        meter.measure([&] {
            auto param = rb.parameter_generator->next(ibody);
            mask.apply(param);
            if (iconstraint == -1) {
                rb.transformer->apply(std::move(param), ibody);
            } else {
                auto constraint = rb.constraints->get_body_constraints(ibody)[iconstraint];
                rb.transformer->apply(std::move(param), constraint);
            }
            rb.transformer->undo(); // keep the molecule in a consistent state between iterations
        });
    };

    // -----------------------------------------------------------------------
    // Hydration shell regeneration.
    // -----------------------------------------------------------------------
    BENCHMARK_ADVANCED("Hydration regeneration") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = make_lar12();
        meter.measure([&] {
            rb.molecule.generate_new_hydration();
        });
    };

    // -----------------------------------------------------------------------
    // Partial histogram update — get_histogram() after a body move.
    // -----------------------------------------------------------------------
    BENCHMARK_ADVANCED("Partial histogram update (get_histogram)") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = make_lar12();
        rb.controller->setup(io::ExistingFile("tests/files/LAR1-2.dat"));
        auto [ibody, iconstraint, mask] = rb.body_selector->next_mask();
        auto param = rb.parameter_generator->next(ibody);
        mask.apply(param);
        if (iconstraint == -1) {
            rb.transformer->apply(std::move(param), ibody);
        } else {
            auto constraint = rb.constraints->get_body_constraints(ibody)[iconstraint];
            rb.transformer->apply(std::move(param), constraint);
        }
        rb.molecule.generate_new_hydration();
        meter.measure([&] {
            return rb.molecule.get_histogram();
        });
    };

    // -----------------------------------------------------------------------
    // Chi2 evaluation (fit_chi2_only) — sinqd cache already warm.
    // -----------------------------------------------------------------------
    BENCHMARK_ADVANCED("fit_chi2_only (sinqd warm)") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = make_lar12();
        rb.controller->setup(io::ExistingFile("tests/files/LAR1-2.dat"));
        auto* fitter = rb.controller->get_fitter();
        fitter->fit_chi2_only(); // prime the sinqd cache
        meter.measure([&] {
            return fitter->fit_chi2_only();
        });
    };

    // -----------------------------------------------------------------------
    // Undo + grid/hydration rebuild (rejected step path).
    // -----------------------------------------------------------------------
    BENCHMARK_ADVANCED("Undo + hydration rebuild") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = make_lar12();
        meter.measure([&] {
            auto [ibody, iconstraint, mask] = rb.body_selector->next_mask();
            auto param = rb.parameter_generator->next(ibody);
            mask.apply(param);
            if (iconstraint == -1) {
                rb.transformer->apply(std::move(param), ibody);
            } else {
                auto constraint = rb.constraints->get_body_constraints(ibody)[iconstraint];
                rb.transformer->apply(std::move(param), constraint);
            }
            rb.transformer->undo();
            rb.molecule.clear_grid();
            rb.molecule.generate_new_hydration();
        });
    };
}

// -----------------------------------------------------------------------
// Symmetry system: SASDJG5 single chain + C2 symmetry (1 body).
// -----------------------------------------------------------------------
TEST_CASE("Rigidbody pipeline step benchmarks: SASDJG5 (C2 symmetry)", "[.][benchmark]") {
    settings::general::verbose = false;
    settings::general::warnings = false;
    settings::molecule::allow_unknown_residues = true;
    settings::molecule::implicit_hydrogens = false;

    BENCHMARK_ADVANCED("Full iteration (prepare + finish)") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = BodySplitter::split("tests/files/SASDJG5_single.pdb");
        rb.molecule.get_body(0).symmetry().add(symmetry::type::p2);
        rb.molecule.generate_new_hydration();
        rb.controller->setup(io::ExistingFile("tests/files/SASDJG5.dat"));
        meter.measure([&] {
            rb.controller->prepare_step();
            rb.controller->finish_step();
        });
    };

    BENCHMARK_ADVANCED("Partial histogram update (get_histogram)") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = BodySplitter::split("tests/files/SASDJG5_single.pdb");
        rb.molecule.get_body(0).symmetry().add(symmetry::type::p2);
        rb.molecule.generate_new_hydration();
        rb.controller->setup(io::ExistingFile("tests/files/SASDJG5.dat"));
        auto [ibody, iconstraint, mask] = rb.body_selector->next_mask();
        auto param = rb.parameter_generator->next(ibody);
        mask.apply(param);
        if (iconstraint == -1) {
            rb.transformer->apply(std::move(param), ibody);
        } else {
            auto constraint = rb.constraints->get_body_constraints(ibody)[iconstraint];
            rb.transformer->apply(std::move(param), constraint);
        }
        rb.molecule.generate_new_hydration();
        meter.measure([&] {
            return rb.molecule.get_histogram();
        });
    };

    BENCHMARK_ADVANCED("fit_chi2_only (sinqd warm)") (Catch::Benchmark::Chronometer meter) {
        Rigidbody rb = BodySplitter::split("tests/files/SASDJG5_single.pdb");
        rb.molecule.get_body(0).symmetry().add(symmetry::type::p2);
        rb.molecule.generate_new_hydration();
        rb.controller->setup(io::ExistingFile("tests/files/SASDJG5.dat"));
        auto* fitter = rb.controller->get_fitter();
        fitter->fit_chi2_only(); // prime the sinqd cache
        meter.measure([&] {
            return fitter->fit_chi2_only();
        });
    };
}
