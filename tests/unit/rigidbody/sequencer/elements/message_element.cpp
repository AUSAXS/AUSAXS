#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/MessageElement.h>
#include <settings/All.h>

#include <iostream>
#include <sstream>
#include <string>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    // the text a print element writes to the console when it runs
    std::string printed(const std::string& message) {
        Sequencer sequencer;
        MessageElement element(&sequencer, message, "white", false);

        bool verbose = settings::general::verbose;
        settings::general::verbose = true;
        std::stringstream out;
        auto* old = std::cout.rdbuf(out.rdbuf());
        element.run();
        std::cout.rdbuf(old);
        settings::general::verbose = verbose;
        return out.str();
    }

    bool contains(const std::string& text, const std::string& part) {return text.find(part) != std::string::npos;}
}

TEST_CASE("MessageElement: placeholders") {
    std::string it = std::to_string(LoopElement::_get_current_iteration());

    CHECK(contains(printed("A: it={iteration} end"), "A: it=" + it + " end"));

    // an unknown placeholder is plain text, and the values after it stay where they were written;
    // the parts used to be paired with the values off by one: "B: {unknown}<it> it= end"
    CHECK(contains(printed("B: {unknown} it={iteration} end"), "B: {unknown} it=" + it + " end"));
    CHECK(contains(printed("D: it={iteration} {x} again={iteration} end"), "D: it=" + it + " {x} again=" + it + " end"));
    CHECK(contains(printed("E: {}{iteration}"), "E: {}" + it));

    // a lone brace is plain text; the text before it used to be printed twice
    std::string lone = printed("C: brace { open");
    CHECK(contains(lone, "C: brace { open"));
    CHECK_FALSE(contains(lone, "C: brace C: brace"));
}

TEST_CASE("MessageElement: an unknown colour is refused when the script is parsed", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 250;
    std::string load =
        "load {\n"
        "    pdb tests/files/2epe.pdb\n"
        "    saxs tests/files/2epe.dat\n"
        "}\n";
    CHECK_NOTHROW(SequenceParser().parse_text(load + "print magenta \"fine\"\n"));
    CHECK_THROWS_AS(SequenceParser().parse_text(load + "print nocolour \"fails at parse time\"\n"), rigidbody::sequencer::except::parse_error);
}
