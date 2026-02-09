#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <data/Molecule.h>
#include <fitter/FitResult.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("Sequencer: programmatic API basic run", "[files]") {
    settings::general::verbose = false;
    settings::grid::min_bins = 500;
    settings::molecule::implicit_hydrogens = false;

    sequencer::Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    auto result = seq
        .setup()
            .load("tests/files/SASDJG5.pdb")
        .end()
        .loop(5)
            .optimize()
        .end()
    .execute();

    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}

TEST_CASE("Sequencer: load with split indices", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    sequencer::Sequencer seq(io::ExistingFile("tests/files/LAR1-2.pdb"));

    // Just verify it doesn't crash during setup with split indices
    REQUIRE_NOTHROW(seq
        .setup()
            .load("tests/files/LAR1-2.pdb", std::vector<int>{9, 99})
        .end()
    );
}

TEST_CASE("Sequencer: load_existing with pre-built Rigidbody", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    auto bodies = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(bodies));

    sequencer::Sequencer seq(io::ExistingFile("tests/files/2epe.dat"));
    auto result = seq
        .setup()
            .load_existing(&rb)
        .end()
        .loop(3)
            .optimize()
        .end()
    .execute();

    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}
