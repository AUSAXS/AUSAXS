#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <settings/All.h>

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
