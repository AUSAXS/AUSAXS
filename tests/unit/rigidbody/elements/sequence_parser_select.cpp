// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/selection/ManualSelect.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/ISymmetry.h>
#include <math/Vector3.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <support/temp_file.h>

#include <vector>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

struct SequenceParserSelectFixture {
    SequenceParserSelectFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        // let the grid size itself to the conformation: the symmetry-isolation case below drives several ungated symmetry perturbations, whose random
        // walk can otherwise carry a symmetry copy outside a grid pinned to this many bins
        settings::grid::min_bins = 25;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config(".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
    }
};

namespace {
    // A symmetry's optimizable parameters, flattened. These are exactly what a ParameterMask governs;
    // a symmetry's other fields (e.g. a cyclic group's rotation angle) are fixed properties of the
    // symmetry type rather than free parameters, and no span exposes them.
    std::vector<double> parameters_of(symmetry::ISymmetry& sym) {
        std::vector<double> flat;
        symmetry::for_each_leaf(sym, [&flat](symmetry::ISymmetry& leaf) {
            for (double t : leaf.span_translation()) {flat.push_back(t);}
            for (double r : leaf.span_rotation()) {flat.push_back(r);}
        });
        return flat;
    }
}

TEST_CASE_METHOD(SequenceParserSelectFixture, "SequenceParser::BodySelectElement inline") {
    SECTION("selecting a body name switches to ManualSelect targeting that body") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select b2\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->body_selector != nullptr);
        CHECK(dynamic_cast<selection::ManualSelect*>(rb->body_selector.get()) != nullptr);
    }

    SECTION("selecting a known strategy name inline works without braces") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select sequential_body\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->body_selector != nullptr);
        CHECK(dynamic_cast<selection::ManualSelect*>(rb->body_selector.get()) == nullptr);
    }

    SECTION("selecting an unknown body name or strategy is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select doesnotexist\n"
        ));
    }

    SECTION("more than one inline argument is rejected") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select b1 b2\n"
        ));
    }

    SECTION("selecting a symmetry tag targets that symmetry alone") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
            "symmetry c3\n"
            "select b1s2\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->body_selector != nullptr);
        REQUIRE(dynamic_cast<selection::ManualSelect*>(rb->body_selector.get()) != nullptr);

        auto selection = rb->body_selector->next_mask();
        CHECK(selection.ibody == 0);
        // a symmetry-targeted selection must never drag its constraint group along, since the group's
        // symmetry parameters are applied to its primary body rather than the selected one
        CHECK(selection.iconstraint == -1);
        CHECK_FALSE(selection.mask.real_translation);
        CHECK_FALSE(selection.mask.real_rotation);
        REQUIRE(selection.mask.target_symmetry.has_value());
        CHECK(selection.mask.target_symmetry.value() == 1); // b1s2 -> the second declared symmetry
    }

    SECTION("a replica tag resolves to the same symmetry as its short form") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
            "symmetry c3\n"
            "select b1s2r2\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->body_selector != nullptr);

        // all replicas of a symmetry are generated from its one shared parameter set, so b1s2r2 targets
        // exactly what b1s2 does
        auto selection = rb->body_selector->next_mask();
        CHECK(selection.ibody == 0);
        REQUIRE(selection.mask.target_symmetry.has_value());
        CHECK(selection.mask.target_symmetry.value() == 1);
    }

    SECTION("a symmetry shared between bodies resolves to the body owning it") {
        // splitting a symmetric body leaves the first fragment owning the shared symmetry and the rest holding non-owning views onto it. Only the owner's
        // copy carries parameters, so a tag naming any other participant has to resolve there - it is the only name that body has for the symmetry.
        auto tag = GENERATE(as<std::string>{}, "b1s1", "b2s1", "b2s1r1");
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/LAR1-2.pdb\n"
            "    saxs tests/files/LAR1-2.dat\n"
            "}\n"
            "symmetry b1 c2\n"
            "split b1 99\n"
            "select " + tag + "\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb->molecule.size_body() == 2);

        auto selection = rb->body_selector->next_mask();
        INFO("selected through tag " << tag);
        CHECK(selection.ibody == 0); // the owning fragment, whichever participant was named
        CHECK(selection.iconstraint == -1);
        REQUIRE(selection.mask.target_symmetry.has_value());
        CHECK(selection.mask.target_symmetry.value() == 0);
    }

    SECTION("an unknown symmetry tag is rejected as a parse error") {
        CHECK_THROWS_AS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
            "select b1s7\n"
        ), sequencer::except::parse_error);

        CHECK_THROWS_AS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
            "select b1s1r99\n"
        ), sequencer::except::parse_error);
    }

    SECTION("named argument block form still works") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "select {\n"
            "    point random_body\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->body_selector != nullptr);
    }

    SECTION("random_symmetry and sequential_symmetry restrict selection to symmetry-owning bodies") {
        // shorthand for `select { body random_body|sequential_body, parameters symmetry }` - the matching body
        // strategy restricted to draws that only ever hit a drivable symmetry, so a parameter_generator with only
        // sym_translate/sym_rotate set never wastes a step on a body that has no symmetry of its own to move
        auto token = GENERATE(as<std::string>{}, "random_symmetry", "sequential_symmetry");
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
            "select " + token + "\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute();
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->body_selector != nullptr);
        CHECK(dynamic_cast<selection::ManualSelect*>(rb->body_selector.get()) == nullptr);

        auto selection = rb->body_selector->next_mask();
        INFO("selected through token " << token);
        CHECK_FALSE(selection.mask.real_translation);
        CHECK_FALSE(selection.mask.real_rotation);
        CHECK(selection.mask.sym_translation);
        CHECK(selection.mask.sym_axis);
        REQUIRE(selection.iconstraint == -1);
    }
}

TEST_CASE_METHOD(SequenceParserSelectFixture, "SequenceParser::BodySelectElement symmetry isolation") {
    auto seq = parse(
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "symmetry c2\n"
        "symmetry c3\n"
        "select b1s2\n"
    );
    REQUIRE(seq != nullptr);
    seq->execute();
    auto rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);
    REQUIRE(rb->body_selector != nullptr);

    auto& body_params = rb->conformation->absolute_parameters.parameters[0];
    REQUIRE(body_params.symmetry_pars.size() == 2);

    auto translation_before = body_params.translation;
    auto rotation_before = body_params.rotation;
    auto untargeted_before = parameters_of(*body_params.symmetry_pars[0]);
    auto targeted_before = parameters_of(*body_params.symmetry_pars[1]);

    // drive several steps exactly as SimpleController::prepare_step does, minus the chi2 accept/reject
    // gate - the gate would revert every rejected step and make the outcome depend on the fit
    for (int i = 0; i < 5; ++i) {
        auto selection = rb->body_selector->next_mask();
        REQUIRE(selection.iconstraint == -1);
        auto par = rb->parameter_generator->next(selection.ibody);
        selection.mask.apply(par);
        rb->transformer->apply(std::move(par), selection.ibody);
    }

    // the targeted symmetry is the only thing that moved: the host body keeps its rigid pose...
    CHECK((body_params.translation - translation_before).magnitude() < 1e-9);
    CHECK((body_params.rotation - rotation_before).magnitude() < 1e-9);

    // ...its other symmetry is frozen...
    CHECK(parameters_of(*body_params.symmetry_pars[0]) == untargeted_before);

    // ...and b1s2 itself is perturbed
    CHECK_FALSE(parameters_of(*body_params.symmetry_pars[1]) == targeted_before);
}

TEST_CASE_METHOD(SequenceParserSelectFixture, "SequenceParser::BodySelectElement symmetry mask split", "[files]") {
    auto mask_of = [this] (const std::string& mask_name) {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
            "select {\n"
            "    body random_body\n"
            "    parameters " + mask_name + "\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        seq->execute(); // the parsed strategy is installed when the element runs
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        return rb->body_selector->next_mask().mask;
    };

    SECTION("symmetry_translation activates the offset alone") {
        auto mask = mask_of("symmetry_translation");
        CHECK_FALSE(mask.real_translation);
        CHECK_FALSE(mask.real_rotation);
        CHECK(mask.sym_translation);
        CHECK_FALSE(mask.sym_axis);
    }

    SECTION("symmetry_axis activates the frame orientation alone") {
        auto mask = mask_of("symmetry_axis");
        CHECK_FALSE(mask.real_translation);
        CHECK_FALSE(mask.real_rotation);
        CHECK_FALSE(mask.sym_translation);
        CHECK(mask.sym_axis);
    }

    SECTION("symmetry still activates both") {
        auto mask = mask_of("symmetry");
        CHECK(mask.sym_translation);
        CHECK(mask.sym_axis);
    }

    SECTION("an unknown mask name is rejected") {
        CHECK_THROWS_AS(mask_of("symmetry_sideways"), sequencer::except::parse_error);
    }
}
