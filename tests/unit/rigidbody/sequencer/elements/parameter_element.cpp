#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <io/ExistingFile.h>
#include <rigidbody/parameters/ParameterAmplitudes.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/ParameterElement.h>
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

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::ParameterElement configuration") {
    SECTION("Configure parameter generation") {
        REQUIRE_NOTHROW(
            parse_sequence(
                "load {\n"
                "    pdb tests/files/2epe.pdb\n"
                "    saxs tests/files/2epe.dat\n"
                "    split 40 80\n"
                "}\n"
                "loop 5\n"
                "    optimize_once\n"
                "end\n"
            )
        );
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

TEST_CASE_METHOD(ParameterParseFixture, "SequenceParser::ParameterElement amplitudes", "[files]") {
    SECTION("translate alone enables only body translation") {
        auto a = amplitudes_of(load() + "parameters {\n    iterations 10\n    translate 5\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0);
        CHECK(a.symmetry_translation == 0);
        CHECK(a.symmetry_rotation == 0);
    }

    SECTION("rotate alone enables only body rotation") {
        auto a = amplitudes_of(load() + "parameters {\n    iterations 10\n    rotate 0.1\n}\n");
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0.1);
        CHECK(a.symmetry_translation == 0);
        CHECK(a.symmetry_rotation == 0);
    }

    SECTION("translate and rotate together no longer drag symmetry along") {
        auto a = amplitudes_of(load_with_symmetry() + "parameters {\n    iterations 10\n    translate 5\n    rotate 0.1\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0.1);
        CHECK(a.symmetry_translation == 0);
        CHECK(a.symmetry_rotation == 0);
    }

    SECTION("sym_translate alone enables only the symmetry offset") {
        auto a = amplitudes_of(load_with_symmetry() + "parameters {\n    iterations 10\n    sym_translate 10\n}\n");
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0);
        CHECK(a.symmetry_translation == 10);
        CHECK(a.symmetry_rotation == 0);
    }

    SECTION("sym_rotate alone enables only the symmetry frame orientation") {
        auto a = amplitudes_of(load_with_symmetry() + "parameters {\n    iterations 10\n    sym_rotate 0.2\n}\n");
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0);
        CHECK(a.symmetry_translation == 0);
        CHECK(a.symmetry_rotation == 0.2);
    }

    SECTION("all four amplitudes can be named at once") {
        auto a = amplitudes_of(load_with_symmetry() +
            "parameters {\n    iterations 10\n    translate 5\n    rotate 0.1\n    sym_translate 10\n    sym_rotate 0.2\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0.1);
        CHECK(a.symmetry_translation == 10);
        CHECK(a.symmetry_rotation == 0.2);
    }

    SECTION("no amplitude at all falls back to the symmetry defaults") {
        auto script = load_with_symmetry() + "parameters {\n    iterations 10\n}\n";
        auto a = amplitudes_of(script);
        CHECK(a.translation == 0);
        CHECK(a.rotation == 0);
        CHECK(a.symmetry_translation > 0);
        CHECK(a.symmetry_rotation == parameter::default_symmetry_rotation());
    }

    SECTION("an explicit zero is the same as leaving the amplitude out") {
        auto a = amplitudes_of(load() + "parameters {\n    iterations 10\n    translate 5\n    rotate 0\n}\n");
        CHECK(a.translation == 5);
        CHECK(a.rotation == 0);
    }

    SECTION("moves are drawn along all axes unless told otherwise") {
        auto a = amplitudes_of(load() + "parameters {\n    iterations 10\n    translate 5\n    rotate 0.1\n}\n");
        CHECK(a.translation_axes == parameter::Axes{});
        CHECK(a.rotation_axes == parameter::Axes{});
    }

    SECTION("axes restrict translations and rotations independently") {
        auto a = amplitudes_of(load() +
            "parameters {\n    iterations 10\n    translate 5\n    rotate 0.1\n    translate_axes xy\n    rotate_axes z\n}\n");
        CHECK(a.translation_axes == parameter::Axes{.x = true, .y = true, .z = false});
        CHECK(a.rotation_axes == parameter::Axes{.x = false, .y = false, .z = true});
    }

    SECTION("the order the axes are written in does not matter") {
        auto a = amplitudes_of(load() + "parameters {\n    iterations 10\n    translate 5\n    translate_axes zx\n}\n");
        CHECK(a.translation_axes == parameter::Axes{.x = true, .y = false, .z = true});
    }
}

TEST_CASE_METHOD(ParameterParseFixture, "SequenceParser::ParameterElement rejects contradictions", "[files]") {
    SECTION("symmetry amplitudes on a molecule without symmetries") {
        CHECK_THROWS_AS(parse(load() + "parameters {\n    iterations 10\n    sym_translate 10\n}\n"), sequencer::except::parse_error);
    }

    SECTION("no amplitude and nothing to fall back on") {
        CHECK_THROWS_AS(parse(load() + "parameters {\n    iterations 10\n}\n"), sequencer::except::parse_error);
    }

    SECTION("mode is no longer an argument") {
        CHECK_THROWS_AS(parse(load() +
            "parameters {\n    iterations 10\n    translate 5\n    mode translate_only\n}\n"), sequencer::except::parse_error);
    }

    SECTION("iterations is still required") {
        CHECK_THROWS_AS(parse(load() + "parameters {\n    translate 5\n}\n"), sequencer::except::parse_error);
    }

    SECTION("axes for a component that is not optimised") {
        CHECK_THROWS_AS(parse(load() + "parameters {\n    iterations 10\n    translate 5\n    rotate_axes z\n}\n"), sequencer::except::parse_error);
        CHECK_THROWS_AS(parse(load() + "parameters {\n    iterations 10\n    rotate 0.1\n    translate_axes xy\n}\n"), sequencer::except::parse_error);
    }

    SECTION("an unknown axis") {
        CHECK_THROWS_AS(parse(load() + "parameters {\n    iterations 10\n    translate 5\n    translate_axes xw\n}\n"), sequencer::except::parse_error);
    }

    SECTION("an axis named twice") {
        CHECK_THROWS_AS(parse(load() + "parameters {\n    iterations 10\n    translate 5\n    translate_axes xx\n}\n"), sequencer::except::parse_error);
    }

}
