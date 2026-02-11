#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <filesystem>
#include <fstream>

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
        
        // Create test output directory
        std::filesystem::create_directories("/tmp/ausaxs_test_output");
    }
    
    ~SequencerElementsFixture() {
        // Cleanup test output directory
        std::filesystem::remove_all("/tmp/ausaxs_test_output");
    }
};

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::SaveElement basic functionality") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Save PDB file") {
        std::string output_path = "/tmp/ausaxs_test_output/test_save.pdb";
        
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(2)
                .optimize()
                .save(output_path)
            .end()
        .execute();
        
        REQUIRE(result != nullptr);
        CHECK(std::filesystem::exists(output_path));
    }
    
    SECTION("Save PDB with counter replacement") {
        std::string output_path = "/tmp/ausaxs_test_output/test_save_%.pdb";
        
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(3)
                .optimize()
                .save(output_path)
            .end()
        .execute();
        
        REQUIRE(result != nullptr);
        // Should create test_save_1.pdb, test_save_2.pdb, test_save_3.pdb
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/test_save_1.pdb"));
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/test_save_2.pdb"));
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/test_save_3.pdb"));
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::EveryNStepElement conditional execution") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Execute every 2 steps") {
        std::string output_path = "/tmp/ausaxs_test_output/every_n_%.pdb";
        
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(5)
                .optimize()
                .every(2)
                    .save(output_path)
                .end()
            .end()
        .execute();
        
        REQUIRE(result != nullptr);
        // Should save on iterations 2, 4
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/every_n_1.pdb"));
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/every_n_2.pdb"));
        CHECK_FALSE(std::filesystem::exists("/tmp/ausaxs_test_output/every_n_3.pdb"));
    }
    
    SECTION("Execute every 3 steps") {
        std::string output_path = "/tmp/ausaxs_test_output/every_3_%.pdb";
        
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(7)
                .optimize()
                .every(3)
                    .save(output_path)
                .end()
            .end()
        .execute();
        
        REQUIRE(result != nullptr);
        // Should save on iterations 3, 6
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/every_3_1.pdb"));
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/every_3_2.pdb"));
        CHECK_FALSE(std::filesystem::exists("/tmp/ausaxs_test_output/every_3_3.pdb"));
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::OnImprovementElement conditional execution") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Execute on improvement only") {
        std::string output_path = "/tmp/ausaxs_test_output/on_improvement_%.pdb";
        
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(5)
                .optimize()
                    .on_improvement()
                        .save(output_path)
                    .end()
                .end()
            .end()
        .execute();
        
        REQUIRE(result != nullptr);
        // At least one improvement should occur
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/on_improvement_1.pdb"));
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::MessageElement") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Simple message without crash") {
        // Just verify it doesn't crash - actual message output is not easily tested
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(2)
                .message("Test message")
                .optimize()
            .end()
        );
    }
    
    SECTION("Message with counter placeholder") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(3)
                .message("Iteration %i")
                .optimize()
            .end()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::OutputFolderElement") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Set output folder") {
        std::string output_folder = "/tmp/ausaxs_test_output/custom_output";
        
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
                .output_folder(output_folder)
            .end()
        );
        
        // Output folder should be created
        CHECK(std::filesystem::exists(output_folder));
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::AutoConstraintsElement") {
    Sequencer seq(io::ExistingFile("tests/files/LAR1-2.pdb"));
    
    SECTION("Generate automatic constraints") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/LAR1-2.pdb", std::vector<int>{9, 99})
                .generate_constraints(
                    settings::rigidbody::ConstraintGenerationStrategyChoice::Linear,
                    settings::rigidbody::DistanceConstraintType::CM,
                    5.0  // overlap strength
                )
            .end()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::RelativeHydrationElement") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Set relative hydration") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
                .relative_hydration(0.9)
            .end()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::ConstraintElement") {
    Sequencer seq(io::ExistingFile("tests/files/LAR1-2.pdb"));
    
    SECTION("Add distance constraint") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/LAR1-2.pdb", std::vector<int>{9, 99})
                .add_constraint(settings::rigidbody::DistanceConstraintType::CM, 0, 1)
            .end()
        );
    }
    
    SECTION("Add overlap constraint") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/LAR1-2.pdb", std::vector<int>{9, 99})
                .add_overlap_constraint(5.0)
            .end()
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::LoopElement nested loops") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Two nested loops") {
        int save_count = 0;
        
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(2)  // outer loop
                .loop(3)  // inner loop
                    .optimize()
                    .save("/tmp/ausaxs_test_output/nested_%.pdb")
                .end()
            .end()
        .execute();
        
        REQUIRE(result != nullptr);
        // Should have 2*3 = 6 saves
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/nested_1.pdb"));
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/nested_6.pdb"));
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
                .parameter(
                    settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly,
                    settings::rigidbody::DecayStrategyChoice::Linear,
                    0.1,  // max rotation
                    0.5,  // max translation
                    0.9   // decay rate
                )
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
                .body_selection(settings::rigidbody::BodySelectStrategyChoice::RandomBodySelect)
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
                .body_selection(settings::rigidbody::BodySelectStrategyChoice::SequentialBodySelect)
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
                .transform(settings::rigidbody::TransformationStrategyChoice::RigidTransform)
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
                .transform(settings::rigidbody::TransformationStrategyChoice::SingleTransform)
                .optimize()
            .end()
        );
    }
}
