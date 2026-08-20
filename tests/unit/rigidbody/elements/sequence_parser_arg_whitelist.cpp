// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/ValidElements.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

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
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config("ausaxs_arg_whitelist_test", ".conf", content);
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

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: unknown named arguments are rejected", "[files]") {
    SECTION("an unknown key in an element that has a whitelist") {
        CHECK_THROWS_AS(parse(load() +
            "parameter {\n"
            "    iterations 10\n"
            "    translate 5\n"
            "    symmetry 1\n"
            "}\n"
        ), sequencer::except::parse_error);
    }

    SECTION("a misspelled key is rejected rather than silently defaulted") {
        CHECK_THROWS_AS(parse(load() +
            "parameter {\n"
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
                "parameter {\n"
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
            "parameter {\n"
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

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: named loops are scoped to one parse", "[files]") {
    auto script = load() +
        "loop L1 10\n"
        "    optimize_once\n"
        "    end\n"
        "end\n";

    CHECK_NOTHROW(parse(script));
    CHECK_NOTHROW(parse(script));
}

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser::SymmetryElement inline forms", "[files]") {
    static const std::string two_bodies =
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n";

    SECTION("bare symmetry name, unambiguous only for a single-body system") {
        auto seq = parse(load() + "symmetry c2\n");
        REQUIRE(seq != nullptr);
        CHECK(seq->_get_rigidbody()->molecule.get_body(0).size_symmetry() == 1);

        CHECK_THROWS_AS(parse(two_bodies + "symmetry c2\n"), sequencer::except::parse_error);
    }

    SECTION("one body and a symmetry") {
        auto seq = parse(two_bodies + "symmetry b2 c2\n");
        REQUIRE(seq != nullptr);
        CHECK(seq->_get_rigidbody()->molecule.get_body(0).size_symmetry() == 0);
        CHECK(seq->_get_rigidbody()->molecule.get_body(1).size_symmetry() == 1);
    }

    SECTION("several bodies share one reference symmetry") {
        auto seq = parse(two_bodies + "symmetry b1 b2 c3\n");
        REQUIRE(seq != nullptr);
        CHECK(seq->_get_rigidbody()->molecule.get_body(0).size_symmetry() == 1);
        CHECK(seq->_get_rigidbody()->molecule.get_body(1).size_symmetry() == 1);
    }

    SECTION("repeated declarations accumulate, replacing the old block form") {
        auto seq = parse(two_bodies +
            "symmetry b1 c2\n"
            "symmetry b2 c3\n"
            "symmetry b1 p2\n"
        );
        REQUIRE(seq != nullptr);
        CHECK(seq->_get_rigidbody()->molecule.get_body(0).size_symmetry() == 2);
        CHECK(seq->_get_rigidbody()->molecule.get_body(1).size_symmetry() == 1);
    }

    SECTION("an unknown body name is rejected") {
        CHECK_THROWS(parse(two_bodies + "symmetry not_a_body c2\n"));
    }

    SECTION("no arguments at all is rejected") {
        CHECK_THROWS_AS(parse(load() + "symmetry\n"), sequencer::except::parse_error);
    }
}

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: inline argument counts are rejected centrally", "[files]") {
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

TEST_CASE_METHOD(ArgWhitelistFixture, "SequenceParser: the two argument forms cannot be combined", "[files]") {
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
