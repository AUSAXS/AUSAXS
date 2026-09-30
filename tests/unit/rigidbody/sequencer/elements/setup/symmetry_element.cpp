#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <io/ExistingFile.h>
#include <math/Vector3.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/parameters/UniformParameterGenerator.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/ValidElements.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/transform/TransformStrategy.h>  // IWYU pragma: keep
#include <settings/All.h>

#include <support/temp_file.h>

#include <algorithm>
#include <string>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;
using namespace ausaxs::rigidbody::constraints;

// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

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

// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

struct SequenceParserSymmetryFixture {
    SequenceParserSymmetryFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
    }

    static std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config(".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
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
        auto* rb = seq->_get_rigidbody();
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
        auto* rb = seq->_get_rigidbody();
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
        auto* rb = seq->_get_rigidbody();
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
        auto* rb = seq->_get_rigidbody();
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
            "symmetry b2 p2-c3\n"
        );
        REQUIRE(seq != nullptr);
        auto* rb = seq->_get_rigidbody();
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
            "symmetry b1 b2 c3\n"
        );
        REQUIRE(seq != nullptr);
        auto* rb = seq->_get_rigidbody();
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

    SECTION("reference symmetry accepts a dihedral base") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry b1 b2 d2\n"
        );
        REQUIRE(seq != nullptr);
        auto* rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);

        REQUIRE(rb->molecule.get_body(0).size_symmetry() == 1);
        REQUIRE(rb->molecule.get_body(1).size_symmetry() == 1);
        auto* ref = dynamic_cast<symmetry::ReferenceSymmetry*>(rb->molecule.get_body(0).symmetry().get(0));
        auto* view = dynamic_cast<symmetry::ReferenceSymmetryView*>(rb->molecule.get_body(1).symmetry().get(0));
        REQUIRE(ref != nullptr);
        REQUIRE(view != nullptr);

        // d2 has order 4 (identity + 3 copies)
        CHECK(ref->repetitions() == 3);
        CHECK(view->repetitions() == 3);
    }

    SECTION("reference symmetry accepts a composite base") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry b1 b2 p2-c3\n"
        );
        REQUIRE(seq != nullptr);
        auto* rb = seq->_get_rigidbody();
        REQUIRE(rb != nullptr);

        REQUIRE(rb->molecule.get_body(0).size_symmetry() == 1);
        auto* ref = dynamic_cast<symmetry::ReferenceSymmetry*>(rb->molecule.get_body(0).symmetry().get(0));
        REQUIRE(ref != nullptr);
        // p2 (inner, 1 copy) nested in c3 (outer, 2 copies) -> (1+1)*(1+2)-1 = 5
        CHECK(ref->repetitions() == 5);

        // for_each_leaf must see through the ReferenceSymmetry into its composite base's two leaves, not treat the wrapper itself as a single (mis-shaped) leaf
        std::vector<symmetry::ISymmetry*> leaves;
        symmetry::for_each_leaf(*ref, [&](symmetry::ISymmetry& leaf) {leaves.push_back(&leaf);});
        REQUIRE(leaves.size() == 2);
    }

    SECTION("symmetry applied to one body does not affect other bodies") {
        auto seq = parse(
            "load {\n"
            "    pdb tests/files/SASDJG5_single.pdb tests/files/SASDJG5_single.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
            "symmetry b2 c2\n"
        );
        REQUIRE(seq != nullptr);
        auto* rb = seq->_get_rigidbody();
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
        "symmetry b1 b2 c3\n"
    );
    REQUIRE(seq != nullptr);
    auto* rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);

    auto* ref = dynamic_cast<symmetry::ReferenceSymmetry*>(rb->molecule.get_body(0).symmetry().get(0));
    auto* view = dynamic_cast<symmetry::ReferenceSymmetryView*>(rb->molecule.get_body(1).symmetry().get(0));
    REQUIRE(ref != nullptr);
    REQUIRE(view != nullptr);

    rigidbody::parameter::UniformParameterGenerator gen(
        rb, 1000, {.symmetry_translation = 5, .symmetry_rotation = 0.5}
    );
    auto nonzero = [](std::span<double> s) {return std::ranges::any_of(s, [](double v) {return v != 0;});};

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
        auto before = view->_get_transform({0, 0, 0}, 1)(probe);

        // transforming the primary body reallocates its symmetry objects; a cached raw pointer would dangle here, but the view re-resolves through the (stable) molecule
        int primary = 0;
        auto params = gen.next(primary);
        rb->transformer->apply(params, primary);

        auto after = view->_get_transform({0, 0, 0}, 1)(probe);
        CHECK((after - before).magnitude() > 1e-6); // the view reflects the updated shared symmetry

        // and it agrees with the primary's current (reallocated) symmetry
        auto* live_ref = dynamic_cast<symmetry::ReferenceSymmetry*>(rb->molecule.get_body(0).symmetry().get(0));
        REQUIRE(live_ref != nullptr);
        auto ref_t = live_ref->_get_transform({0, 0, 0}, 1)(probe);
        CHECK((after - ref_t).magnitude() < 1e-9);
    }
}
