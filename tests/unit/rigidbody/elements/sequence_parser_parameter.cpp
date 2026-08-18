// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/ParameterElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/parameters/ParameterAmplitudes.h>
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

struct ParameterParseFixture {
    ParameterParseFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config("ausaxs_parameter_parse_test", ".conf", content);
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
    parameter::ParameterAmplitudes amplitudes_of(const std::string& script) {
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

TEST_CASE_METHOD(ParameterParseFixture, "SequenceParser::ParameterElement amplitudes", "[files]") {
    SECTION("translate alone enables only body translation") {
        auto a = amplitudes_of(load() + "parameter {\n    iterations 10\n    translate 5\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0);
        CHECK(a.symmetry_translation == 0);
        CHECK(a.symmetry_rotation == 0);
    }

    SECTION("rotate alone enables only body rotation") {
        auto a = amplitudes_of(load() + "parameter {\n    iterations 10\n    rotate 0.1\n}\n");
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0.1);
        CHECK(a.symmetry_translation == 0);
        CHECK(a.symmetry_rotation == 0);
    }

    SECTION("translate and rotate together no longer drag symmetry along") {
        auto a = amplitudes_of(load_with_symmetry() + "parameter {\n    iterations 10\n    translate 5\n    rotate 0.1\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0.1);
        CHECK(a.symmetry_translation == 0);
        CHECK(a.symmetry_rotation == 0);
    }

    SECTION("sym_translate alone enables only the symmetry offset") {
        auto a = amplitudes_of(load_with_symmetry() + "parameter {\n    iterations 10\n    sym_translate 10\n}\n");
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0);
        CHECK(a.symmetry_translation == 10);
        CHECK(a.symmetry_rotation == 0);
    }

    SECTION("sym_rotate alone enables only the symmetry frame orientation") {
        auto a = amplitudes_of(load_with_symmetry() + "parameter {\n    iterations 10\n    sym_rotate 0.2\n}\n");
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0);
        CHECK(a.symmetry_translation == 0);
        CHECK(a.symmetry_rotation == 0.2);
    }

    SECTION("all four amplitudes can be named at once") {
        auto a = amplitudes_of(load_with_symmetry() +
            "parameter {\n    iterations 10\n    translate 5\n    rotate 0.1\n    sym_translate 10\n    sym_rotate 0.2\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0.1);
        CHECK(a.symmetry_translation == 10);
        CHECK(a.symmetry_rotation == 0.2);
    }

    SECTION("no amplitude at all falls back to the symmetry defaults") {
        auto script = load_with_symmetry() + "parameter {\n    iterations 10\n}\n";
        auto a = amplitudes_of(script);
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0);
        CHECK(a.symmetry_translation > 0);
        CHECK(a.symmetry_rotation == parameter::default_symmetry_rotation());
    }

    SECTION("an explicit zero is the same as leaving the amplitude out") {
        auto a = amplitudes_of(load() + "parameter {\n    iterations 10\n    translate 5\n    rotate 0\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0);
    }
}

TEST_CASE_METHOD(ParameterParseFixture, "SequenceParser::ParameterElement modes", "[files]") {
    SECTION("mode symmetry retargets translate and rotate onto the symmetry components") {
        auto a = amplitudes_of(load_with_symmetry() +
            "parameter {\n    iterations 10\n    mode symmetry\n    translate 10\n    rotate 0.2\n}\n");
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0);
        CHECK(a.symmetry_translation == 10);
        CHECK(a.symmetry_rotation == 0.2);
    }

    SECTION("mode symmetry with only translate gives symmetry translation alone") {
        auto a = amplitudes_of(load_with_symmetry() + "parameter {\n    iterations 10\n    mode symmetry\n    translate 10\n}\n");
        CHECK(a.symmetry_translation == 10);
        CHECK(a.symmetry_rotation == 0);
    }

    SECTION("mode symmetry with only rotate gives symmetry rotation alone") {
        auto a = amplitudes_of(load_with_symmetry() + "parameter {\n    iterations 10\n    mode symmetry\n    rotate 0.2\n}\n");
        CHECK(a.symmetry_translation == 0);
        CHECK(a.symmetry_rotation == 0.2);
    }

    SECTION("mode symmetry without amplitudes uses the defaults") {
        auto a = amplitudes_of(load_with_symmetry() + "parameter {\n    iterations 10\n    mode symmetry\n}\n");
        CHECK(a.symmetry_translation > 0);
        CHECK(a.symmetry_rotation == parameter::default_symmetry_rotation());
    }

    SECTION("symmetry_only remains accepted as a spelling of symmetry") {
        auto a = amplitudes_of(load_with_symmetry() + "parameter {\n    iterations 10\n    mode symmetry_only\n    translate 10\n}\n");
        CHECK(a.symmetry_translation == 10);
    }

    SECTION("translate_only") {
        auto a = amplitudes_of(load() + "parameter {\n    iterations 10\n    mode translate_only\n    translate 5\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0);
    }

    SECTION("rotate_only") {
        auto a = amplitudes_of(load() + "parameter {\n    iterations 10\n    mode rotate_only\n    rotate 0.1\n}\n");
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0.1);
    }

    SECTION("both") {
        auto a = amplitudes_of(load() + "parameter {\n    iterations 10\n    mode both\n    translate 5\n    rotate 0.1\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0.1);
    }

    SECTION("strategy is an alias for mode") {
        auto a = amplitudes_of(load() + "parameter {\n    iterations 10\n    strategy rotate_only\n    rotate 0.1\n}\n");
        CHECK(a.rotation == 0.1);
    }
}

TEST_CASE_METHOD(ParameterParseFixture, "SequenceParser::ParameterElement rejects contradictions", "[files]") {
    SECTION("mode symmetry alongside an explicit symmetry amplitude") {
        CHECK_THROWS_AS(parse(load_with_symmetry() +
            "parameter {\n    iterations 10\n    mode symmetry\n    sym_translate 10\n}\n"), sequencer::except::parse_error);
    }

    SECTION("a mode rejects the amplitudes it has no use for") {
        CHECK_THROWS_AS(parse(load() +
            "parameter {\n    iterations 10\n    mode translate_only\n    translate 5\n    rotate 0.1\n}\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load_with_symmetry() +
            "parameter {\n    iterations 10\n    mode both\n    translate 5\n    rotate 0.1\n    sym_rotate 0.2\n}\n"), sequencer::except::parse_error);
    }

    SECTION("a mode still requires its own amplitude") {
        CHECK_THROWS_AS(parse(load() + "parameter {\n    iterations 10\n    mode translate_only\n}\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "parameter {\n    iterations 10\n    mode rotate_only\n}\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "parameter {\n    iterations 10\n    mode both\n    translate 5\n}\n"), sequencer::except::parse_error);
    }

    SECTION("symmetry amplitudes on a molecule without symmetries") {
        CHECK_THROWS_AS(parse(load() + "parameter {\n    iterations 10\n    sym_translate 10\n}\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "parameter {\n    iterations 10\n    mode symmetry\n}\n"), sequencer::except::parse_error);
    }

    SECTION("no amplitude and nothing to fall back on") {
        CHECK_THROWS_AS(parse(load() + "parameter {\n    iterations 10\n}\n"), sequencer::except::parse_error);
    }

    SECTION("iterations is still required") {
        CHECK_THROWS_AS(parse(load() + "parameter {\n    translate 5\n}\n"), sequencer::except::parse_error);
    }

    SECTION("an unknown mode") {
        CHECK_THROWS_AS(parse(load() + "parameter {\n    iterations 10\n    mode sideways\n    translate 5\n}\n"), sequencer::except::parse_error);
    }
}
