#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/elements/DebugCounterElement.h>
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
    // collects every debug counter in the element tree, in script order
    void find_counters(LoopElement* loop, std::vector<DebugCounterElement*>& out) {
        for (auto& e : loop->_get_elements()) {
            if (auto* counter = dynamic_cast<DebugCounterElement*>(e.get())) {out.push_back(counter);}
            else if (auto* nested = dynamic_cast<LoopElement*>(e.get())) {find_counters(nested, out);}
        }
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::EveryNStepElement conditional execution") {
    auto steps = GENERATE(1, 2, 3, 5);
    auto seq = SequenceParser().parse_text(
        "load {\n"
        "    pdb tests/files/2epe.pdb\n"
        "    saxs tests/files/2epe.dat\n"
        "    split 40 80\n"
        "}\n"
        "loop 6\n"
        "    debug_counter outer\n"
        "    every " + std::to_string(steps) + "\n"
        "        debug_counter inner\n"
        "    end\n"
        "end\n"
    );
    REQUIRE(seq != nullptr);
    REQUIRE_NOTHROW(seq->execute());

    std::vector<DebugCounterElement*> counters;
    find_counters(seq.get(), counters);
    REQUIRE(counters.size() == 2);
    CHECK(counters[0]->name() == "outer");
    CHECK(counters[0]->count() == 6);
    CHECK(counters[1]->name() == "inner");
    CHECK(counters[1]->count() == 6/steps);
}
