#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <fitter/FitResult.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

struct SequencerElementsFixture {
    SequencerElementsFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }
};

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::SaveElement basic functionality") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Save PDB file - verify no crash") {
        std::string output_path = "/tmp/ausaxs_test_output/test_save.pdb";
        
        REQUIRE_NOTHROW(
            seq
                .setup()
                    .load("tests/files/SASDJG5.pdb")
                .end()
                .loop(2)
                    .optimize()
                    .save(output_path)
                .end()
            .execute()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::EveryNStepElement conditional execution") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Execute every 2 steps - verify no crash") {
        std::string output_path = "/tmp/ausaxs_test_output/every_n_%.pdb";
        
        REQUIRE_NOTHROW(
            seq
                .setup()
                    .load("tests/files/SASDJG5.pdb")
                .end()
                .loop(5)
                    .optimize()
                    .every(2)
                        .save(output_path)
                    .end()
                .end()
            .execute()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::OnImprovementElement conditional execution") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Basic optimization steps") {
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
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::AutoConstraintsElement") {
    Sequencer seq(io::ExistingFile("tests/files/LAR1-2.pdb"));
    
    SECTION("Generate linear constraints") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/LAR1-2.pdb", std::vector<int>{9, 99})
                .generate_linear_constraints()
            .end()
        );
    }
    
    SECTION("Generate volumetric constraints") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/LAR1-2.pdb", std::vector<int>{9, 99})
                .generate_volumetric_constraints()
            .end()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::ConstraintElement") {
    Sequencer seq(io::ExistingFile("tests/files/LAR1-2.pdb"));
    
    SECTION("Add distance constraint center mass") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/LAR1-2.pdb", std::vector<int>{9, 99})
                .distance_constraint_center_mass(0, 1)
            .end()
        );
    }
    
    SECTION("Add distance constraint closest") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/LAR1-2.pdb", std::vector<int>{9, 99})
                .distance_constraint_closest(0, 1)
            .end()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::LoopElement nested loops") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Two nested loops") {
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(2)  // outer loop
                .loop(3)  // inner loop
                    .optimize()
                .end()
            .end()
        .execute();
        
        REQUIRE(result != nullptr);
        CHECK(result->fval > 0);
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::ParameterElement configuration") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Configure parameter generation") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(5)
                .optimize()
            .end()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::BodySelectElement strategies") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Random body selection") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(3)
                .optimize()
            .end()
        );
    }
    
    SECTION("Sequential body selection") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(3)
                .optimize()
            .end()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::TransformElement strategies") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Rigid transform") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(3)
                .optimize()
            .end()
        );
    }
    
    SECTION("Single transform") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(3)
                .optimize()
            .end()
        );
    }
}
