// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <string>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

struct LoopCountFixture {
    LoopCountFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    static std::string load() {
        return
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n";
    }

    // the total step count the given script body would report, as Sequencer::execute calculates it
    unsigned int total_iterations_of(const std::string& body) {
        test::TempFile config("ausaxs_loop_count_test", ".conf", load() + body);
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
