// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <io/ExistingFile.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/ValidElements.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <string>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

struct ArgWhitelistFixture {
    ArgWhitelistFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
    }

    static std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config(".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
    }

    // every script needs a loaded molecule before any other element can be parsed
    static std::string load() {
        return
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n";
    }
};

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: unknown named arguments are rejected") {
    SECTION("an unknown key in an element that has a whitelist") {
        CHECK_THROWS_AS(parse(load() +
            "parameters {\n"
            "    iterations 10\n"
            "    translate 5\n"
            "    symmetry 1\n"
            "}\n"
        ), sequencer::except::parse_error);
    }

    SECTION("a misspelled key is rejected rather than silently defaulted") {
        CHECK_THROWS_AS(parse(load() +
            "parameters {\n"
            "    iteration 10\n"
            "    translate 5\n"
            "}\n"
        ), sequencer::except::parse_error);
    }

    SECTION("an unknown key in an element that takes no named arguments") {
        CHECK_THROWS_AS(parse(load() +
            "select {\n"
            "    nonsense random_body\n"
            "}\n"
        ), sequencer::except::parse_error);
    }

    SECTION("the error names the offending key and lists the valid ones") {
        try {
            parse(load() +
                "parameters {\n"
                "    iterations 10\n"
                "    translate 5\n"
                "    symmetry 1\n"
                "}\n"
            );
            FAIL("expected a parse_error");
        } catch (const sequencer::except::parse_error& e) {
            std::string what = e.what();
            CHECK(what.find("\"symmetry\"") != std::string::npos);
            CHECK(what.find("\"iterations\"") != std::string::npos); // the valid-argument list
        }
    }

    SECTION("every valid key is still accepted") {
        CHECK_NOTHROW(parse(load() + "symmetry c2\n" +
            "parameters {\n"
            "    iterations 10\n"
            "    translate 5\n"
            "    rotate 0.1\n"
            "    sym_translate 10\n"
            "    sym_rotate 0.2\n"
            "    decay exponential\n"
            "}\n"
        ));
    }

    SECTION("aliases of valid keys are accepted") {
        CHECK_NOTHROW(parse(load() +
            "select {\n"
            "    body random_body\n"
            "    mask symmetry\n" // alias of "parameters"
            "}\n"
        ));
    }
}

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: named loops are scoped to one parse") {
    auto script = load() +
        "loop L1 10\n"
        "    optimize_once\n"
        "    end\n"
        "end\n";

    CHECK_NOTHROW(parse(script));
    CHECK_NOTHROW(parse(script));

    SECTION("a duplicate name within one script says so") {
        try {
            parse(script + "loop L1 2\n    optimize_once\n    end\nend\n");
            FAIL("expected a parse_error");
        } catch (const sequencer::except::parse_error& e) {
            CHECK(std::string(e.what()).find("already exists") != std::string::npos);
        }
    }
}

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: inline argument counts are rejected centrally") {
    SECTION("too few inline arguments") {
        CHECK_THROWS_AS(parse(load() + "rename b1\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "merge b1\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "delete\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "symmetry\n"), sequencer::except::parse_error);
    }

    SECTION("too many inline arguments") {
        CHECK_THROWS_AS(parse(load() + "save a.pdb b.pdb\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "loop 2\nend surplus\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "update structure extra\n"), sequencer::except::parse_error);
    }

    SECTION("an element taking no inline arguments rejects them") {
        // optimize_once/optimize_step takes no inline arguments at all, unlike overlap_strength (which gained a
        // 1-2 argument inline form and so no longer fits this case)
        CHECK_THROWS_AS(parse(load() + "optimize_once 5 10\n"), sequencer::except::parse_error);
    }

    SECTION("the error names the expected arguments") {
        try {
            parse(load() + "rename b1\n");
            FAIL("expected a parse_error");
        } catch (const sequencer::except::parse_error& e) {
            std::string what = e.what();
            CHECK(what.find("[old name] [new name]") != std::string::npos);
            CHECK(what.find("but got 1") != std::string::npos);
        }
    }

    SECTION("open-ended elements accept long argument lists") {
        CHECK_NOTHROW(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "merge b1 b2 b3\n"
        ));
    }
}

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: the two argument forms cannot be combined") {
    // tokens between the element name and the brace used to be discarded without a word, so this was previously accepted
    SECTION("an inline argument in front of a block") {
        CHECK_THROWS_AS(parse(load() + "print red {\n    msg \"hello\"\n}\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "output somewhere {\n    mode relative\n}\n"), sequencer::except::parse_error);
    }

    SECTION("the error says the forms cannot be combined") {
        try {
            parse(load() + "print red {\n    msg \"hello\"\n}\n");
            FAIL("expected a parse_error");
        } catch (const sequencer::except::parse_error& e) {
            CHECK(std::string(e.what()).find("Cannot combine inline and named arguments") != std::string::npos);
        }
    }

    SECTION("either form on its own is still fine") {
        CHECK_NOTHROW(parse(load() + "print red \"hello\"\n"));
        CHECK_NOTHROW(parse(load() + "print {\n    msg \"hello\"\n    colour red\n}\n"));
    }
}

TEST_CASE("SequenceParser: every element declares a well-formed inline signature") {
    using namespace ausaxs::rigidbody::sequencer::detail;
    for (int i = 0; i < static_cast<int>(ElementType::COUNT); ++i) {
        auto type = static_cast<ElementType>(i);
        auto signature = valid_inline_arguments(type);
        INFO("element type index " << i);
        CHECK(signature.min <= signature.max);
        CHECK(signature.max <= unbounded_inline_args);
        for (const auto& name : signature.names) {
            CHECK_FALSE(name.empty());
        }
        // a name for every slot the element can actually be given, so error messages and the GUI never show a blank
        if (signature.max == 0) {
            CHECK(signature.names.empty());
        } else {
            CHECK_FALSE(signature.names.empty());
        }
    }
}

TEST_CASE("SequenceParser: every element declares its own named arguments") {
    // guards against an element gaining a named argument in _parse without adding it to _valid_arguments, which the central
    // whitelist would then reject. Only checks that the declaration exists and is well-formed.
    using namespace ausaxs::rigidbody::sequencer::detail;
    for (int i = 0; i < static_cast<int>(ElementType::COUNT); ++i) {
        auto type = static_cast<ElementType>(i);
        auto args = valid_arguments(type);
        INFO("element type index " << i);
        for (const auto& arg : args) {
            CHECK_FALSE(arg.empty());
        }
    }
}

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: an unknown element reports its own line") {
    // the load block is lines 1-4; the unknown element opens a block on line 5 which ends on line 7
    try {
        parse(load() + "parameter {\n    iterations 10\n}\n");
        FAIL("expected a parse_error");
    } catch (const sequencer::except::parse_error& e) {
        CHECK(std::string(e.what()).find("(line 5)") != std::string::npos);
        CHECK(std::string(e.what()).find("Did you mean \"parameters\"?") != std::string::npos);
    }
}

TEST_CASE("SequenceParser: element names must match a registered name exactly") {
    using namespace ausaxs::rigidbody::sequencer::detail;
    CHECK(get_type("loop") == ElementType::LoopBegin);
    CHECK(get_type("parameters") == ElementType::Parameter);

    // a registered name followed by anything else is a different, unknown name: "loop3" used to parse as a bare "loop"
    CHECK_THROWS_AS(get_type("loop3"), sequencer::except::parse_error);
    CHECK_THROWS_AS(get_type("loop_3"), sequencer::except::parse_error);
    CHECK_THROWS_AS(get_type("save-all"), sequencer::except::parse_error);
    CHECK_THROWS_AS(get_type("parameters_strategy"), sequencer::except::parse_error);

    // the message names the line and the closest valid name
    auto message = [] (std::string_view name, int line) -> std::string {
        try {get_type(name, line);} catch (const sequencer::except::parse_error& e) {return e.what();}
        return "";
    };
    CHECK(message("parameter", 7).find("(line 7)") != std::string::npos);
    CHECK(message("parameter", 7).find("Did you mean \"parameters\"?") != std::string::npos);
    CHECK(message("optimise_once", -1).find("Did you mean \"optimize_once\"?") != std::string::npos);
    CHECK(message("parameter", -1).find("line") == std::string::npos);
    CHECK(message("xyzzy", -1).find("Did you mean") == std::string::npos);

    // every name the parser accepts is advertised to the GUI
    for (const auto& name : valid_elements()) {
        INFO(name);
        CHECK_NOTHROW(get_type(name));
    }
}

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: overlap_strength sets the overlap penalty") {
    struct OverlapConstraint : constraints::OverlapConstraint {using constraints::OverlapConstraint::weight;}; // weight() is protected
    using Catch::Matchers::WithinRel;
    const double default_weight = std::exp(-5*0.5);

    SECTION("a scaling alone multiplies the default weight") {
        auto seq = parse(load() + "overlap_strength 2\n");
        CHECK_THAT(OverlapConstraint::weight(0.5), WithinRel(2*default_weight, 1e-12));
    }

    SECTION("a scaling and a distance give a quadratic well of that range") {
        auto seq = parse(load() + "overlap_strength 2 5\n");
        CHECK_THAT(OverlapConstraint::weight(2.5), WithinRel(0.5, 1e-12));
        CHECK(OverlapConstraint::weight(5) == 0);
    }

    SECTION("the block form agrees with the inline form") {
        auto seq = parse(load() + "overlap_strength {\n    scaling 2\n    distance 5\n}\n");
        CHECK_THAT(OverlapConstraint::weight(2.5), WithinRel(0.5, 1e-12));
    }

    SECTION("a distance alone uses unit scaling") {
        auto seq = parse(load() + "overlap_strength {\n    distance 5\n}\n");
        CHECK_THAT(OverlapConstraint::weight(2.5), WithinRel(0.25, 1e-12));
    }

    SECTION("malformed values are rejected") {
        CHECK_THROWS_AS(parse(load() + "overlap_strength\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "overlap_strength x\n"), sequencer::except::parse_error);
        CHECK_THROWS(parse(load() + "overlap_strength nan\n"));
        CHECK_THROWS(parse(load() + "overlap_strength 2 5abc\n"));
    }
}

TEST_CASE("SequenceParser: script numbers are parsed strictly") {
    using search::parse_int, search::parse_double;
    CHECK(parse_int("12") == 12);
    CHECK(parse_int("+12") == 12);
    CHECK(parse_int("-3") == -3);
    CHECK(parse_double("2.5") == 2.5);
    CHECK(parse_double("1e3") == 1000);
    CHECK(parse_double("-.5") == -0.5);

    // std::stoi read these as 1, 2 and 3, and std::stod accepted nan and inf
    for (const auto* bad : {"1e3", "2.5", "3abc", "", "+", "99999999999999999999"}) {
        INFO(bad);
        CHECK_THROWS(parse_int(bad));
    }
    for (const auto* bad : {"nan", "inf", "-inf", "+nan", "1e999", "1.5x", ""}) {
        INFO(bad);
        CHECK_THROWS(parse_double(bad));
    }
}

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: numeric arguments are checked") {
    auto parameter = [] (const std::string& body) {return "parameters {\n" + body + "}\n";};
    std::string P = parameter("    iterations 4\n    translate 1\n    rotate 1\n");

    SECTION("loop counts are positive integers") {
        auto loop_iterations = [&] (const std::string& line) {
            auto seq = parse(load() + P + line + "\n    optimize_once\n    end\nend\n");
            auto* loop = dynamic_cast<LoopElement*>(seq->_get_elements().back().get());
            REQUIRE(loop != nullptr);
            return loop->_get_loop_iterations();
        };
        CHECK(loop_iterations("loop 3") == 3);
        CHECK(loop_iterations("loop") == 4);
        CHECK(loop_iterations("loop outer") == 4);
        CHECK(loop_iterations("loop outer 7") == 7);

        // these ran 1, 2, 3, 0 and 0 steps
        for (std::string bad : {"loop 1e3", "loop 2.5", "loop 3abc", "loop -3", "loop 0", "loop outer 2.5", "loop outer -1"}) {
            INFO(bad);
            CHECK_THROWS_AS(parse(load() + P + bad + "\n    optimize_once\n    end\nend\n"), sequencer::except::parse_error);
        }
    }

    SECTION("a duplicate loop name is reported as such") {
        try {
            parse(load() + P + "loop a 2\n    optimize_once\n    end\nend\nloop a 2\n    optimize_once\n    end\nend\n");
            FAIL("expected a parse_error");
        } catch (const sequencer::except::parse_error& e) {
            CHECK(std::string(e.what()).find("already exists") != std::string::npos);
        }
    }

    SECTION("every takes a positive step count") {
        CHECK_NOTHROW(parse(load() + P + "loop 4\n    every 2\n        optimize_once\n        end\n    end\nend\n"));
        // every 0 was a SIGFPE in the step count and in run()
        for (std::string bad : {"every 0", "every -1", "every 1.5"}) {
            INFO(bad);
            CHECK_THROWS(parse(load() + P + "loop 4\n    " + bad + "\n        optimize_once\n        end\n    end\nend\n"));
        }
    }

    SECTION("parameter iterations and amplitudes") {
        for (std::string bad : {
            "    iterations 0\n    translate 1\n",
            "    iterations -2\n    translate 1\n",
            "    iterations 1e3\n    translate 1\n",
            "    iterations 4\n    translate nan\n",
            "    iterations 4\n    rotate inf\n",
            "    iterations 4\n    translate -1\n"
        }) {
            INFO(bad);
            CHECK_THROWS(parse(load() + parameter(bad)));
        }
    }

    SECTION("seeds and split indices are whole integers") {
        CHECK_THROWS_AS(parse("seed 1e3\n" + load()), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "split b1 9.5\n"), sequencer::except::parse_error);
    }
}
