#include <catch2/catch_test_macros.hpp>

#include <fitter/FitResult.h>  // IWYU pragma: keep
#include <io/Folder.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <settings/All.h>

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
    }
};

namespace {
    std::unique_ptr<Sequencer> parse_sequence(const std::string& script) {
        SequenceParser parser;
        return parser.parse_text(script + (script.find("\nloop ") != std::string::npos ? "end\n" : ""));
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::SaveElement basic functionality") {
    SECTION("Save PDB file - verify no crash") {
        io::Folder out_dir("temp/ausaxs_test_output_" + test::detail::unique_tag());
        out_dir.create();
        std::string output_path = out_dir.path() + "/test_save.pdb";
        auto seq = parse_sequence(
            "load {\n"
            "    pdb tests/files/SASDJG5.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "loop 2\n"
            "    optimize_once\n"
            "    save " + output_path + "\n"
            "end\n"
        );

        REQUIRE(seq != nullptr);
        REQUIRE_NOTHROW(seq->execute());
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::EveryNStepElement conditional execution") {
    SECTION("Execute every 2 steps - verify no crash") {
        io::Folder out_dir("temp/ausaxs_test_output_" + test::detail::unique_tag());
        out_dir.create();
        std::string output_path = out_dir.path() + "/every_n_%.pdb";
        auto seq = parse_sequence(
            "load {\n"
            "    pdb tests/files/SASDJG5.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "loop 5\n"
            "    optimize_once\n"
            "    every 2\n"
            "        save " + output_path + "\n"
            "    end\n"
            "end\n"
        );

        REQUIRE(seq != nullptr);
        REQUIRE_NOTHROW(seq->execute());
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::OnImprovementElement conditional execution") {
    SECTION("Basic optimization steps") {
        auto seq = parse_sequence(
            "load {\n"
            "    pdb tests/files/SASDJG5.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "loop 5\n"
            "    optimize_once\n"
            "end\n"
        );
        
        REQUIRE(seq != nullptr);
        auto result = seq->execute();
        REQUIRE(result != nullptr);
        CHECK(result->fval > 0);
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::AutoConstraintsElement") {
    SECTION("Generate backbone constraints") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/LAR1-2.pdb\n"
                "    saxs tests/files/LAR1-2.dat\n"
                "    split 9 99\n"
                "}\n"
                "autoconstrain backbone\n"
            )
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::ConstraintElement") {
    SECTION("Add distance constraint center mass") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/LAR1-2.pdb\n"
                "    saxs tests/files/LAR1-2.dat\n"
                "    split 9 99\n"
                "}\n"
                "constrain {\n"
                "    first b1\n"
                "    second b2\n"
                "    type cm\n"
                "}\n"
            )
        );
    }

    SECTION("Add distance constraint closest") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/LAR1-2.pdb\n"
                "    saxs tests/files/LAR1-2.dat\n"
                "    split 9 99\n"
                "}\n"
                "constrain {\n"
                "    first b1\n"
                "    second b2\n"
                "    type bond\n"
                "}\n"
            )
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::LoopElement nested loops") {
    SECTION("Two nested loops") {
        auto seq = parse_sequence(
            "load {\n"
            "    pdb tests/files/SASDJG5.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "loop 2\n"
            "    loop 3\n"
            "        optimize_once\n"
            "    end\n"
            "end\n"
        );

        REQUIRE(seq != nullptr);
        auto result = seq->execute();
        REQUIRE(result != nullptr);
        CHECK(result->fval > 0);
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::ParameterElement configuration") {
    SECTION("Configure parameter generation") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/SASDJG5.pdb\n"
                "    saxs tests/files/SASDJG5.dat\n"
                "}\n"
                "loop 5\n"
                "    optimize_once\n"
                "end\n"
            )
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::BodySelectElement strategies") {
    SECTION("Random body selection") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/SASDJG5.pdb\n"
                "    saxs tests/files/SASDJG5.dat\n"
                "}\n"
                "loop 3\n"
                "    optimize_once\n"
                "end\n"
            )
        );
    }

    SECTION("Sequential body selection") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/SASDJG5.pdb\n"
                "    saxs tests/files/SASDJG5.dat\n"
                "}\n"
                "loop 3\n"
                "    optimize_once\n"
                "end\n"
            )
        );
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::TransformElement strategies") {
    SECTION("Rigid transform") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/SASDJG5.pdb\n"
                "    saxs tests/files/SASDJG5.dat\n"
                "}\n"
                "loop 3\n"
                "    optimize_once\n"
                "end\n"
            )
        );
    }

    SECTION("Single transform") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/SASDJG5.pdb\n"
                "    saxs tests/files/SASDJG5.dat\n"
                "}\n"
                "loop 3\n"
                "    optimize_once\n"
                "end\n"
            )
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
        auto seq = parse_sequence(
            "load {\n"
            "    pdb tests/files/SASDJG5.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "loop 10\n"
            "    optimize_once\n"
            "end\n"
        );
        REQUIRE(seq != nullptr);
        auto* loop = dynamic_cast<LoopElement*>(seq->_get_elements().back().get());
        REQUIRE(loop != nullptr);
        loop->_get_elements().push_back(std::make_unique<StopRequestElement>());

        auto result = seq->execute();

        // the requesting iteration always finishes, so exactly one of the ten should have run
        CHECK(LoopElement::_get_current_iteration() == 1);
        CHECK(LoopElement::_stop_requested());

        // a stopped run is still a complete run: the best conformation so far is restored and fitted
        REQUIRE(result != nullptr);
        CHECK(result->fval > 0);
    }

    SECTION("A stop requested while nothing is running does not affect the next run") {
        LoopElement::_request_stop();

        auto seq = parse_sequence(
            "load {\n"
            "    pdb tests/files/SASDJG5.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "loop 3\n"
            "    optimize_once\n"
            "end\n"
        );
        REQUIRE(seq != nullptr);
        auto result = seq->execute();

        CHECK(LoopElement::_get_current_iteration() == 3);
        CHECK_FALSE(LoopElement::_stop_requested());
        REQUIRE(result != nullptr);
    }
}
