#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fitter/FitResult.h>  // IWYU pragma: keep
#include <io/ExistingFile.h>
#include <io/Folder.h>
#include <rigidbody/parameters/ParameterAmplitudes.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/elements/ParameterElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <string>

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

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::LoopElement nested loops") {
    SECTION("Two nested loops") {
        auto seq = parse_sequence(
            "load {\n"
            "    pdb tests/files/2epe.pdb\n"
            "    saxs tests/files/2epe.dat\n"
            "    split 40 80\n"
            "}\n"
            "loop 2\n"
            "    loop 2\n"
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
            "    pdb tests/files/2epe.pdb\n"
            "    saxs tests/files/2epe.dat\n"
            "    split 40 80\n"
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
            "    pdb tests/files/2epe.pdb\n"
            "    saxs tests/files/2epe.dat\n"
            "    split 40 80\n"
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

// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

struct ParameterParseFixture {
    ParameterParseFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
    }

    static std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config(".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
    }

    static std::string load() {
        return
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n";
    }

    static std::string load_with_symmetry() {return load() + "symmetry c2\n";}

    // the amplitudes the parsed "parameter" element ended up with
    static parameter::ParameterAmplitudes amplitudes_of(const std::string& script) {
        auto seq = parse(script);
        REQUIRE(seq != nullptr);
        for (auto& element : seq->_get_elements()) {
            if (auto* p = dynamic_cast<ParameterElement*>(element.get())) {
                return p->get_parameter_strategy()->get_amplitudes();
            }
        }
        FAIL("script contained no parameter element");
        return {};
    }
};

TEST_CASE_METHOD(ParameterParseFixture, "SequenceParser::LoopElement iteration deduction", "[files]") {
    // a bare "loop" deduces its iteration count from the last parameter element, walking up the owner chain to find one.
    // With no parameter element anywhere, that walk reaches the Sequencer, which has no owner to continue to.
    SECTION("a bare loop with no parameter element to deduce from is a parse error") {
        CHECK_THROWS_AS(parse(load() + "loop\n    optimize_once\n    end\nend\n"), sequencer::except::parse_error);
    }

    SECTION("a named bare loop fails the same way") {
        CHECK_THROWS_AS(parse(load() + "loop outer\n    optimize_once\n    end\nend\n"), sequencer::except::parse_error);
    }

    SECTION("a nested bare loop fails the same way") {
        CHECK_THROWS_AS(
            parse(load() + "loop 5\n    loop\n        optimize_once\n        end\n    end\nend\n"),
            sequencer::except::parse_error
        );
    }

    SECTION("an explicit iteration count needs no parameter element") {
        CHECK_NOTHROW(parse(load() + "loop 5\n    optimize_once\n    end\nend\n"));
    }

    SECTION("a parameter element supplies the count to a bare loop") {
        auto seq = parse(load() + "parameter {\n    iterations 7\n    translate 5\n}\nloop\n    optimize_once\n    end\nend\n");
        REQUIRE(seq != nullptr);
        for (auto& element : seq->_get_elements()) {
            if (auto* loop = dynamic_cast<LoopElement*>(element.get())) {
                CHECK(loop->_get_loop_iterations() == 7);
                return;
            }
        }
        FAIL("script contained no loop element");
    }
}

// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

struct LoopCountFixture {
    LoopCountFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
    }

    static std::string load() {
        return
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n";
    }

    // the total step count the given script body would report, as Sequencer::execute calculates it
    static int total_iterations_of(const std::string& body) {
        test::TempFile config(".conf", load() + body);
        SequenceParser parser;
        auto seq = parser.parse_file(config);
        REQUIRE(seq != nullptr);
        LoopElement::_recount_total_iterations(seq.get());
        return LoopElement::_get_total_iterations();
    }
};

TEST_CASE_METHOD(LoopCountFixture, "SequenceParser::LoopElement total iteration count", "[files]") {
    SECTION("a single loop counts its own iterations") {
        CHECK(total_iterations_of(
            "loop 5\n"
            "    optimize_once\n"
            "    end\n"
            "end\n"
        ) == 5);
    }

    SECTION("a loop containing only another loop does not count itself") {
        CHECK(total_iterations_of(
            "loop 5\n"
            "    loop 50\n"
            "        optimize_once\n"
            "        end\n"
            "    end\n"
            "end\n"
        ) == 250);
    }

    SECTION("a loop can both optimize itself and contain a nested loop") {
        CHECK(total_iterations_of(
            "loop 5\n"
            "    optimize_once\n"
            "    end\n"
            "    loop 10\n"
            "        optimize_once\n"
            "        end\n"
            "    end\n"
            "end\n"
        ) == 55);
    }

    SECTION("sibling loops are summed before being multiplied by their parent") {
        CHECK(total_iterations_of(
            "loop 5\n"
            "    loop 50\n"
            "        optimize_once\n"
            "        end\n"
            "    end\n"
            "    loop 20\n"
            "        optimize_once\n"
            "        end\n"
            "    end\n"
            "end\n"
        ) == 350);
    }

    SECTION("a copy loop counts as much as its target") {
        CHECK(total_iterations_of(
            "loop 5\n"
            "    loop L1 50\n"
            "        optimize_once\n"
            "        end\n"
            "    end\n"
            "    loop copy L1\n"
            "end\n"
        ) == 500);
    }

    SECTION("steps inside an every-n block are only run every nth iteration") {
        CHECK(total_iterations_of(
            "loop 10\n"
            "    every 2\n"
            "        optimize_once\n"
            "        end\n"
            "    end\n"
            "end\n"
        ) == 5);
    }

    SECTION("a loop with no optimization steps at all counts nothing") {
        CHECK(total_iterations_of(
            "loop 5\n"
            "    save dummy.pdb\n"
            "end\n"
        ) == 0);
    }
}
