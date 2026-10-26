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
    CHECK(contains(printed("B: it={iteration} again={iteration} end"), "B: it=" + it + " again=" + it + " end"));
}

TEST_CASE("MessageElement: an unknown placeholder or stray brace is refused") {
    // the user evidently meant to print a value there, so it must not silently become text
    Sequencer sequencer;
    CHECK_THROWS_AS(MessageElement(&sequencer, "{unknown}", "white", false), rigidbody::sequencer::except::parse_error);
    CHECK_THROWS_AS(MessageElement(&sequencer, "it={iteration} {x}", "white", false), rigidbody::sequencer::except::parse_error);
    CHECK_THROWS_AS(MessageElement(&sequencer, "{}", "white", false), rigidbody::sequencer::except::parse_error);
    CHECK_THROWS_AS(MessageElement(&sequencer, "{chi}", true), rigidbody::sequencer::except::parse_error);

    // braces are reserved for placeholders, so a stray one is refused too
    CHECK_THROWS_AS(MessageElement(&sequencer, "brace { open", "white", false), rigidbody::sequencer::except::parse_error);
    CHECK_THROWS_AS(MessageElement(&sequencer, "brace } close", "white", false), rigidbody::sequencer::except::parse_error);
    CHECK_THROWS_AS(MessageElement(&sequencer, "{ {iteration}", "white", false), rigidbody::sequencer::except::parse_error);
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
