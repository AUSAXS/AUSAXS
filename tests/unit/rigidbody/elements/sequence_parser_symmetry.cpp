// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/parameters/ParameterGenerationStrategies.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <math/Vector3.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <algorithm>
#include <fstream>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;
using namespace ausaxs::rigidbody::constraints;

struct SequenceParserSymmetryFixture {
    SequenceParserSymmetryFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        static int counter = 0;
        std::string path = "/tmp/ausaxs_seq_sym_test_" + std::to_string(counter++) + ".conf";
        std::ofstream f(path);
        f << content;
        f.close();
        SequenceParser parser;
        return parser.parse_file(path);
    }
};

TEST_CASE_METHOD(SequenceParserSymmetryFixture, "SequenceParser::SymmetryElement") {
    SECTION("c2 symmetry creates one symmetry on the body") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 1);
    }

    SECTION("no symmetry directive leaves body without symmetries") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 0);
    }

    SECTION("two symmetry directives produce two symmetries on the body") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry c2\n"
            "symmetry c2\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 2);
    }

    SECTION("composite symmetry p2-c3 builds a nested CompositeSymmetry") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry p2-c3\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        REQUIRE(rb->molecule.get_body(0).size_symmetry() == 1);

        auto* comp = dynamic_cast<symmetry::CompositeSymmetry*>(rb->molecule.get_body(0).symmetry().get(0));
        REQUIRE(comp != nullptr);
        // p2 (inner, 1 copy) nested in c3 (outer, 2 copies) -> (1+1)*(1+2)-1 = 5
        CHECK(comp->repetitions() == 5);
    }

    SECTION("composite symmetry can be applied to a named body in a block") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry {\n"
            "    b2 p2-c3\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 0);
        REQUIRE(rb->molecule.get_body(1).size_symmetry() == 1);
        auto* comp = dynamic_cast<symmetry::CompositeSymmetry*>(rb->molecule.get_body(1).symmetry().get(0));
        REQUIRE(comp != nullptr);
        // p2 (inner, 1 copy) nested in c3 (outer, 2 copies) -> (1+1)*(1+2)-1 = 5
        CHECK(comp->repetitions() == 5);
    }

    SECTION("reference symmetry shares one symmetry across several bodies") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry {\n"
            "    bodies \"b1 b2\" c3\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);

        // the primary body owns a ReferenceSymmetry; the other holds a non-owning view of it
        REQUIRE(rb->molecule.get_body(0).size_symmetry() == 1);
        REQUIRE(rb->molecule.get_body(1).size_symmetry() == 1);
        auto* ref = dynamic_cast<symmetry::ReferenceSymmetry*>(rb->molecule.get_body(0).symmetry().get(0));
        auto* view = dynamic_cast<symmetry::ReferenceSymmetryView*>(rb->molecule.get_body(1).symmetry().get(0));
        REQUIRE(ref != nullptr);
        REQUIRE(view != nullptr);

        // the view forwards to the primary's symmetry, so both report the same repetitions (c3 -> 2)
        CHECK(ref->repetitions() == 2);
        CHECK(view->repetitions() == 2);
        CHECK(view->target() == ref);
    }

    SECTION("reference symmetry rejects a non-cyclic base") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry {\n"
            "    bodies \"b1 b2\" t\n"
            "}\n"
        ));
    }

    SECTION("symmetry applied to one body does not affect other bodies") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry {\n"
            "    b2 c2\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        auto rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);
        CHECK(rb->molecule.get_body(0).size_symmetry() == 0);
        CHECK(rb->molecule.get_body(1).size_symmetry() == 1);
    }
}

TEST_CASE_METHOD(SequenceParserSymmetryFixture, "SequenceParser: reference symmetry refinement") {
    auto seq = parse(
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "symmetry {\n"
        "    bodies \"b1 b2\" c3\n"
        "}\n"
    );
    REQUIRE(seq != nullptr);
    auto rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);

    auto* ref = dynamic_cast<symmetry::ReferenceSymmetry*>(rb->molecule.get_body(0).symmetry().get(0));
    auto* view = dynamic_cast<symmetry::ReferenceSymmetryView*>(rb->molecule.get_body(1).symmetry().get(0));
    REQUIRE(ref != nullptr);
    REQUIRE(view != nullptr);

    rigidbody::parameter::SymmetryOnly gen(rb, settings::rigidbody::iterations, 5, 0.5);
    auto nonzero = [](std::span<double> s) {return std::any_of(s.begin(), s.end(), [](double v) {return v != 0;});};

    SECTION("the shared symmetry is optimisable, the view is inert") {
        // the primary body's symmetry is perturbed...
        auto p_primary = gen.next(0);
        REQUIRE(p_primary.symmetry_pars.has_value());
        REQUIRE(p_primary.symmetry_pars.value().size() == 1);
        auto* delta = dynamic_cast<symmetry::ReferenceSymmetry*>(p_primary.symmetry_pars.value()[0].get());
        REQUIRE(delta != nullptr);
        CHECK((nonzero(delta->span_translation()) || nonzero(delta->span_rotation())));

        // ...but the view contributes no optimisable parameters of its own
        auto p_view = gen.next(1);
        REQUIRE(p_view.symmetry_pars.has_value());
        REQUIRE(p_view.symmetry_pars.value().size() == 1);
        auto* view_delta = dynamic_cast<symmetry::ReferenceSymmetryView*>(p_view.symmetry_pars.value()[0].get());
        REQUIRE(view_delta != nullptr);
        CHECK(view_delta->span_translation().empty());
        CHECK(view_delta->span_rotation().empty());
    }

    SECTION("a view survives transformation of the primary body and tracks it") {
        Vector3<double> probe{1, 2, 3};
        auto before = view->get_transform({0, 0, 0}, 1)(probe);

        // transforming the primary body reallocates its symmetry objects; a cached raw pointer
        // would dangle here, but the view re-resolves through the (stable) molecule
        unsigned int primary = 0;
        auto params = gen.next(primary);
        rb->transformer->apply(std::move(params), primary);

        auto after = view->get_transform({0, 0, 0}, 1)(probe);
        CHECK((after - before).magnitude() > 1e-6); // the view reflects the updated shared symmetry

        // and it agrees with the primary's current (reallocated) symmetry
        auto* live_ref = dynamic_cast<symmetry::ReferenceSymmetry*>(rb->molecule.get_body(0).symmetry().get(0));
        REQUIRE(live_ref != nullptr);
        auto ref_t = live_ref->get_transform({0, 0, 0}, 1)(probe);
        CHECK((after - ref_t).magnitude() < 1e-9);
    }
}

TEST_CASE_METHOD(SequenceParserSymmetryFixture, "SequenceParser::ConstraintElement real-real") {
    auto seq = parse(
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "constrain {\n"
        "    body1 b1\n"
        "    body2 b2\n"
        "    type attract\n"
        "    distance 30\n"
        "}\n"
    );
    REQUIRE(seq != nullptr);
    auto rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);
    // non_discoverable_constraints[0] is the pre-added OverlapConstraint; ours is at the back
    REQUIRE(rb->constraints->non_discoverable_constraints.size() >= 2);
    auto* c = dynamic_cast<IDistanceConstraint*>(rb->constraints->non_discoverable_constraints.back().get());
    REQUIRE(c != nullptr);

    SECTION("ibody1 and ibody2 reference the two distinct bodies") {
        CHECK(c->ibody1 == 0);
        CHECK(c->ibody2 == 1);
    }

    SECTION("both isym values indicate real (non-symmetry) bodies") {
        CHECK(c->isym1 == std::make_pair(-1, 0));
        CHECK(c->isym2 == std::make_pair(-1, 0));
    }
}

TEST_CASE_METHOD(SequenceParserSymmetryFixture, "SequenceParser::ConstraintElement real-symmetry") {
    auto seq = parse(
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "symmetry c2\n"
        "constrain {\n"
        "    body1 b1s1\n"
        "    body2 b1\n"
        "    type attract\n"
        "    distance 30\n"
        "}\n"
    );
    REQUIRE(seq != nullptr);
    auto rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);
    REQUIRE(rb->constraints->non_discoverable_constraints.size() >= 2);
    auto* c = dynamic_cast<IDistanceConstraint*>(rb->constraints->non_discoverable_constraints.back().get());
    REQUIRE(c != nullptr);

    SECTION("both constrained atoms belong to body 0") {
        CHECK(c->ibody1 == 0);
        CHECK(c->ibody2 == 0);
    }

    SECTION("body1 isym tracks the first symmetry, first replica") {
        CHECK(c->isym1 == std::make_pair(0, 1));
    }

    SECTION("body2 isym indicates a real (non-symmetry) body") {
        CHECK(c->isym2 == std::make_pair(-1, 0));
    }
}

TEST_CASE_METHOD(SequenceParserSymmetryFixture, "SequenceParser::ConstraintElement symmetry-symmetry") {
    auto seq = parse(
        "load {\n"
        "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "symmetry {\n"
        "    b1 c2\n"
        "    b2 c2\n"
        "}\n"
        "constrain {\n"
        "    body1 b1s1\n"
        "    body2 b2s1\n"
        "    type attract\n"
        "    distance 50\n"
        "}\n"
    );
    REQUIRE(seq != nullptr);
    auto rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);
    REQUIRE(rb->constraints->non_discoverable_constraints.size() >= 2);
    auto* c = dynamic_cast<IDistanceConstraint*>(rb->constraints->non_discoverable_constraints.back().get());
    REQUIRE(c != nullptr);

    SECTION("body indices reference the two distinct bodies") {
        CHECK(c->ibody1 == 0);
        CHECK(c->ibody2 == 1);
    }

    SECTION("both isym values track the first symmetry, first replica of their respective bodies") {
        CHECK(c->isym1 == std::make_pair(0, 1));
        CHECK(c->isym2 == std::make_pair(0, 1));
    }
}
