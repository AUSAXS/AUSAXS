// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <support/temp_file.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

struct SequenceParserQuotingFixture {
    SequenceParserQuotingFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
        settings::general::output = "output/";
    }

    ~SequenceParserQuotingFixture() {
        settings::general::output = "output/";
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config(".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
    }

    static constexpr const char* load_block =
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
    ;
};

// A quoted path is the only way to write a path containing spaces, and on Windows such a path
// commonly ends in a '\' separator. The closing quote must therefore still terminate the token,
// or the separator-plus-quote is absorbed into the path and every later write fails ('"' is not
// a legal character in a Windows file name).
TEST_CASE_METHOD(SequenceParserQuotingFixture, "SequenceParser: quoted paths") {
    SECTION("quoted path with spaces") {
        auto seq = parse(std::string(load_block) + "output \"out dir/sub dir/\"\n");
        REQUIRE(seq != nullptr);
        CHECK(settings::general::output == "out dir/sub dir/");
    }

    SECTION("quoted windows path ending in a separator") {
        auto seq = parse(std::string(load_block) + "output \"C:\\out dir\\sub dir\\\"\n");
        REQUIRE(seq != nullptr);
        CHECK(settings::general::output == "C:\\out dir\\sub dir\\");
    }

    SECTION("quoted windows path without a trailing separator gains one") {
        auto seq = parse(std::string(load_block) + "output \"C:\\out dir\\sub dir\"\n");
        REQUIRE(seq != nullptr);
        CHECK(settings::general::output == "C:\\out dir\\sub dir/");
    }

    SECTION("trailing comment after a quoted windows path") {
        auto seq = parse(std::string(load_block) + "output \"C:\\out dir\\\" # where the results go\n");
        REQUIRE(seq != nullptr);
        CHECK(settings::general::output == "C:\\out dir\\");
    }

    SECTION("unquoted path is unaffected") {
        auto seq = parse(std::string(load_block) + "output outdir/\n");
        REQUIRE(seq != nullptr);
        CHECK(settings::general::output == "outdir/");
    }
}
