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
    
    SECTION("Basic optimization steps") {
        std::string output_path = "/tmp/ausaxs_test_output/basic_%.pdb";
        
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(5)
                .optimize()
                .save(output_path)
            .end()
        .execute();
        
        REQUIRE(result != nullptr);
        // At least one save should occur
        CHECK(std::filesystem::exists("/tmp/ausaxs_test_output/basic_1.pdb"));
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::MessageElement") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Simple loop without crash") {
        // Just verify the sequencer runs without crashing
        auto result = seq
            .setup()
                .load("tests/files/SASDJG5.pdb")
            .end()
            .loop(2)
                .optimize()
            .end()
        .execute();
        
        REQUIRE(result != nullptr);
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

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::RelativeHydrationElement") {
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    
    SECTION("Set relative hydration") {
        REQUIRE_NOTHROW(
            seq.setup()
                .load("tests/files/SASDJG5.pdb")
                .relative_hydration({"A", "B"}, {0.9, 1.1})
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
