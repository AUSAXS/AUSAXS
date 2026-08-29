#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <fitter/FitResult.h>
#include <settings/All.h>
#include <io/ExistingFile.h>
#include <io/Folder.h>

#include <support/temp_file.h>

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
        io::Folder out_dir("temp/ausaxs_test_output_" + test::detail::unique_tag());
        out_dir.create();
        std::string output_path = out_dir.path() + "/test_save.pdb";

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
        io::Folder out_dir("temp/ausaxs_test_output_" + test::detail::unique_tag());
        out_dir.create();
        std::string output_path = out_dir.path() + "/every_n_%.pdb";

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
    
    SECTION("Generate backbone constraints") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/LAR1-2.pdb", std::vector<int>{9, 99})
                .generate_backbone_constraints()
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

namespace {
    // requests a stop the first time it is run, so the surrounding loop should not start another iteration
    struct StopRequestElement : GenericElement {
        void run() override {LoopElement::_request_stop();}
    };
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::LoopElement stop request") {
    SECTION("Stop request ends the loop after the current iteration") {
        Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
        auto& loop = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(10);
        loop.optimize();
        loop._get_elements().push_back(std::make_unique<StopRequestElement>());

        auto result = loop.end().execute();

        // the requesting iteration always finishes, so exactly one of the ten should have run
        CHECK(LoopElement::_get_current_iteration() == 1);
        CHECK(LoopElement::_stop_requested());

        // a stopped run is still a complete run: the best conformation so far is restored and fitted
        REQUIRE(result != nullptr);
        CHECK(result->fval > 0);
    }

    SECTION("A stop requested while nothing is running does not affect the next run") {
        LoopElement::_request_stop();

        Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(3)
                .optimize()
            .end()
        .execute();

        CHECK(LoopElement::_get_current_iteration() == 3);
        CHECK_FALSE(LoopElement::_stop_requested());
        REQUIRE(result != nullptr);
    }
}
