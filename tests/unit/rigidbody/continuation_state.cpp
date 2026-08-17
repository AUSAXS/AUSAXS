// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <rigidbody/continuation/ContinuationState.h>
#include <rigidbody/constraints/AttractorConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/constraints/RepellerConstraint.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/atoms/AtomMetadata.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <settings/All.h>

#include <support/temp_file.h>

#include <fstream>
#include <memory>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::continuation;

namespace {
    // A scratch path for a state file. TempFile only exists to own the name and delete it again; the contents are
    // written by write_continuation_state itself.
    struct StateFile {
        test::TempFile file{"continuation", std::string(continuation_extension)};
        operator const io::File&() const {return file;}
    };

    std::unique_ptr<Rigidbody> make_rigidbody(std::vector<Body> bodies) {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
        return std::make_unique<Rigidbody>(Molecule{std::move(bodies)});
    }

    // Three bodies of two carbons each, spaced far enough apart not to overlap.
    std::vector<Body> three_bodies() {
        std::vector<Body> bodies;
        for (int i = 0; i < 3; ++i) {
            bodies.emplace_back(std::vector{
                AtomFF({10.0*i,       0, 0}, form_factor::form_factor_t::C),
                AtomFF({10.0*i + 1.5, 0, 0}, form_factor::form_factor_t::C)
            });
        }
        return bodies;
    }
}

TEST_CASE("ContinuationState::atoms_round_trip") {
    auto rb = make_rigidbody(three_bodies());
    StateFile path;

    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    REQUIRE(state.molecule.size_body() == rb->molecule.size_body());
    for (unsigned int i = 0; i < rb->molecule.size_body(); ++i) {
        const auto& expected = rb->molecule.get_body(i).get_atoms();
        const auto& actual = state.molecule.get_body(i).get_atoms();
        REQUIRE(actual.size() == expected.size());
        for (unsigned int j = 0; j < expected.size(); ++j) {
            // exact equality, not approximate: a continuation must resume from precisely the pose it stored, so any
            // loss here would silently perturb the system between runs
            CHECK(actual[j].coordinates() == expected[j].coordinates());
            CHECK(actual[j].weight() == expected[j].weight());
            CHECK(actual[j].form_factor_type() == expected[j].form_factor_type());
        }
    }
}

TEST_CASE("ContinuationState::hydration_is_not_carried_over") {
    // Hydration is derived data, not rigid-body state: every refinement calls generate_new_hydration(), which clears
    // the layer and rebuilds it from the grid. Carrying the previous shell across does not survive that, but it is
    // present while the grid is first built, which perturbs the shell that gets regenerated - enough to shift chi2 for
    // a structure that has not moved at all. So it is deliberately dropped.
    std::vector<Body> bodies;
    bodies.emplace_back(
        std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)},
        std::vector{Water({1, 2, 3}), Water({4, 5, 6})}
    );
    auto rb = make_rigidbody(std::move(bodies));
    REQUIRE(rb->molecule.get_body(0).size_water() == 2);

    StateFile path;
    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    CHECK(state.molecule.get_body(0).size_water() == 0);
    CHECK(state.molecule.get_body(0).get_atoms().size() == 1);
}

TEST_CASE("ContinuationState::metadata_round_trip") {
    // The whole reason this format exists rather than reusing .pdb: writing a molecule to .pdb discards the residue
    // numbering and backbone flags, leaving a continued run unable to re-split by residue or generate backbone bonds.
    std::vector<Body> bodies;
    bodies.emplace_back(std::vector{
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        AtomFF({1, 0, 0}, form_factor::form_factor_t::N),
        AtomFF({2, 0, 0}, form_factor::form_factor_t::C)
    });
    AtomMetadata md;
    md.backbone = std::vector{backbone_t::c_alpha, backbone_t::none, backbone_t::c_alpha};
    md.residue_seq = std::vector{7, 7, 8};
    md.occupancy = std::vector{1.0f, 0.5f, 0.25f};
    bodies[0].set_metadata(std::move(md));

    auto rb = make_rigidbody(std::move(bodies));
    StateFile path;

    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    const auto& restored = state.molecule.get_body(0).get_metadata();
    REQUIRE(restored.has_value());
    REQUIRE(restored->backbone.has_value());
    REQUIRE(restored->residue_seq.has_value());
    REQUIRE(restored->occupancy.has_value());
    CHECK(*restored->backbone == std::vector{backbone_t::c_alpha, backbone_t::none, backbone_t::c_alpha});
    CHECK(*restored->residue_seq == std::vector{7, 7, 8});
    CHECK(*restored->occupancy == std::vector{1.0f, 0.5f, 0.25f});
}

TEST_CASE("ContinuationState::absent_metadata_stays_absent") {
    auto rb = make_rigidbody(three_bodies());
    StateFile path;

    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    CHECK_FALSE(state.molecule.get_body(0).get_metadata().has_value());
}

TEST_CASE("ContinuationState::body_names_round_trip") {
    auto rb = make_rigidbody(three_bodies());
    StateFile path;

    write_continuation_state(path, *rb, {"first", "b2", "renamed"});
    auto state = read_continuation_state(path);

    CHECK(state.body_names == std::vector<std::string>{"first", "b2", "renamed"});
}

TEST_CASE("ContinuationState::symmetry_parameters_round_trip") {
    auto rb = make_rigidbody(three_bodies());
    rb->molecule.get_body(0).symmetry().add(symmetry::create("c4"));

    // perturb the parameters, standing in for what a refinement does to them: a symmetry restored at its *default*
    // pose would throw away everything the previous run optimised
    auto translation = rb->molecule.get_body(0).symmetry().get(0)->span_translation();
    for (std::size_t i = 0; i < translation.size(); ++i) {translation[i] = 1.5 + i;}
    auto rotation = rb->molecule.get_body(0).symmetry().get(0)->span_rotation();
    for (std::size_t i = 0; i < rotation.size(); ++i) {rotation[i] = 0.25 * (i + 1);}

    StateFile path;
    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    REQUIRE(state.molecule.get_body(0).symmetry().get().size() == 1);
    REQUIRE(state.molecule.get_body(1).symmetry().get().size() == 0);
    auto restored = state.molecule.get_body(0).symmetry().get(0);
    CHECK(restored->type_name() == "c4");
    CHECK(restored->repetitions() == rb->molecule.get_body(0).symmetry().get(0)->repetitions());

    auto restored_translation = restored->span_translation();
    REQUIRE(restored_translation.size() == translation.size());
    for (std::size_t i = 0; i < translation.size(); ++i) {CHECK(restored_translation[i] == translation[i]);}

    auto restored_rotation = restored->span_rotation();
    REQUIRE(restored_rotation.size() == rotation.size());
    for (std::size_t i = 0; i < rotation.size(); ++i) {CHECK(restored_rotation[i] == rotation[i]);}
}

TEST_CASE("ContinuationState::composite_symmetry_round_trips") {
    auto rb = make_rigidbody(three_bodies());
    rb->molecule.get_body(0).symmetry().add(
        std::make_unique<symmetry::CompositeSymmetry>(symmetry::create("p2"), symmetry::create("c3"))
    );
    auto expected_name = rb->molecule.get_body(0).symmetry().get(0)->type_name();

    StateFile path;
    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    auto restored = state.molecule.get_body(0).symmetry().get(0);
    CHECK(restored->type_name() == expected_name);
    CHECK(dynamic_cast<const symmetry::CompositeSymmetry*>(restored) != nullptr);
}

TEST_CASE("ContinuationState::stacked_symmetries_round_trip") {
    auto rb = make_rigidbody(three_bodies());
    rb->molecule.get_body(1).symmetry().add(symmetry::create("c2"));
    rb->molecule.get_body(1).symmetry().add(symmetry::create("c3"));

    StateFile path;
    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    REQUIRE(state.molecule.get_body(1).symmetry().get().size() == 2);
    CHECK(state.molecule.get_body(1).symmetry().get(0)->type_name() == "c2");
    CHECK(state.molecule.get_body(1).symmetry().get(1)->type_name() == "c3");
}

TEST_CASE("ContinuationState::shared_symmetry_round_trips") {
    // A symmetry shared across bodies is one owned ReferenceSymmetry plus a view on every other participant. Restoring
    // the participants as independent symmetries would quietly turn one fitted parameter set into several.
    auto rb = make_rigidbody(three_bodies());
    rb->molecule.get_body(0).symmetry().add(
        std::make_unique<symmetry::ReferenceSymmetry>(symmetry::create("c2"), std::vector<int>{0, 1}, std::vector<int>{0, 0}, &rb->molecule)
    );
    rb->molecule.get_body(1).symmetry().add(std::make_unique<symmetry::ReferenceSymmetryView>(&rb->molecule, 0, 0));

    StateFile path;
    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    auto* owner = dynamic_cast<const symmetry::ReferenceSymmetry*>(state.molecule.get_body(0).symmetry().get(0));
    REQUIRE(owner != nullptr);
    CHECK(owner->bodies == std::vector<int>{0, 1});
    CHECK(owner->slots == std::vector<int>{0, 0});
    CHECK(owner->base->type_name() == "c2");

    auto* view = dynamic_cast<const symmetry::ReferenceSymmetryView*>(state.molecule.get_body(1).symmetry().get(0));
    REQUIRE(view != nullptr);
    CHECK(view->primary_body == 0);
    CHECK(view->symmetry_index == 0);
}

TEST_CASE("ContinuationState::constraints_round_trip") {
    auto rb = make_rigidbody(three_bodies());
    using Kind = ContinuationConstraint::Kind;

    // built through the restoring constructors so the stored values are exactly what the test asserts on, rather than
    // whatever the deriving constructors would compute from this synthetic geometry
    rb->constraints->add_constraint(std::make_unique<constraints::DistanceConstraintBond>(
        constraints::restore, &rb->molecule, 0, 1, 1, 0, std::pair{-1, -1}, std::pair{-1, -1}, 3.5));
    rb->constraints->add_constraint(std::make_unique<constraints::DistanceConstraintCM>(
        constraints::restore, &rb->molecule, 0, 0, 2, 1, std::pair{-1, -1}, std::pair{-1, -1}, 21.0));
    rb->constraints->add_constraint(std::make_unique<constraints::AttractorConstraint>(
        constraints::restore, &rb->molecule, 1, 0, 2, 0, std::pair{0, 1}, std::pair{-1, -1}, 30.0));
    rb->constraints->add_constraint(std::make_unique<constraints::RepellerConstraint>(
        constraints::restore, &rb->molecule, 0, 1, 2, 1, std::pair{-1, -1}, std::pair{1, 2}, 45.0));
    rb->constraints->add_constraint(std::make_unique<constraints::DistanceConstraintAtom>(
        constraints::restore, &rb->molecule, 1, 1, 2, 0, std::pair{-1, -1}, std::pair{-1, -1}, 12.5));

    StateFile path;
    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    REQUIRE(state.constraints.size() == 5);

    CHECK(state.constraints[0].kind == Kind::bond);
    CHECK(state.constraints[0].ibody1 == 0);
    CHECK(state.constraints[0].iatom1 == 1);
    CHECK(state.constraints[0].ibody2 == 1);
    CHECK(state.constraints[0].iatom2 == 0);
    CHECK(state.constraints[0].d_target == 3.5);

    CHECK(state.constraints[1].kind == Kind::cm);
    CHECK(state.constraints[1].d_target == 21.0);

    CHECK(state.constraints[2].kind == Kind::attractor);
    CHECK(state.constraints[2].isym1 == std::pair{0, 1});
    CHECK(state.constraints[2].d_target == 30.0);

    CHECK(state.constraints[3].kind == Kind::repeller);
    CHECK(state.constraints[3].isym2 == std::pair{1, 2});
    CHECK(state.constraints[3].d_target == 45.0);

    CHECK(state.constraints[4].kind == Kind::atom);
    CHECK(state.constraints[4].d_target == 12.5);
}

TEST_CASE("ContinuationState::stores_non_discoverable_constraints") {
    // ConstraintManager files everything except backbone bonds under non_discoverable_constraints - that split is about
    // per-body lookup during refinement, not about what a constraint is. Reading only the discoverable list would drop
    // every cm/attract/repel constraint a user declared, so this pins both lists being written.
    auto rb = make_rigidbody(three_bodies());
    rb->constraints->add_constraint(std::make_unique<constraints::DistanceConstraintCM>(
        constraints::restore, &rb->molecule, 0, 0, 1, 0, std::pair{-1, -1}, std::pair{-1, -1}, 10.0));
    REQUIRE(rb->constraints->discoverable_constraints.empty());

    StateFile path;
    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    REQUIRE(state.constraints.size() == 1);
    CHECK(state.constraints[0].kind == ContinuationConstraint::Kind::cm);
    CHECK(state.constraints[0].d_target == 10.0);
}

TEST_CASE("ContinuationState::skips_constraints_that_are_not_distance_constraints") {
    // A Rigidbody creates its own overlap constraint, which is not a distance constraint and is regenerated by whoever
    // loads the state. Writing it out would double it up on every continuation.
    auto rb = make_rigidbody(three_bodies());
    StateFile path;

    write_continuation_state(path, *rb, {});
    auto state = read_continuation_state(path);

    CHECK(state.constraints.empty());
}

TEST_CASE("ContinuationState::restored_constraints_keep_their_target") {
    // The deriving constructors set d_target from the *current* geometry. Re-deriving on load would therefore reset
    // every target to the distance the bodies had already drifted to, silently neutralising the constraint.
    auto rb = make_rigidbody(three_bodies());
    constraints::DistanceConstraintCM derived(&rb->molecule, 0, 1);
    constraints::DistanceConstraintCM restored(
        constraints::restore, &rb->molecule, 0, derived.iatom1, 1, derived.iatom2, {-1, -1}, {-1, -1}, 4.0);

    CHECK(derived.d_target != 4.0);
    CHECK(restored.d_target == 4.0);
    CHECK(restored.evaluate() != derived.evaluate());
}

TEST_CASE("ContinuationState::restored_bond_ignores_the_distance_limit") {
    // A refinement may stretch a backbone bond well past the 4 Å a freshly declared bond is allowed to span. A state
    // that was legal to write must stay legal to read, so the restoring path must not re-apply that check.
    auto rb = make_rigidbody(three_bodies());
    CHECK_NOTHROW(constraints::DistanceConstraintBond(
        constraints::restore, &rb->molecule, 0, 0, 2, 0, std::pair{-1, -1}, std::pair{-1, -1}, 20.0));
}

TEST_CASE("ContinuationState::rejects_a_foreign_file") {
    test::TempFile foreign("not_a_state", std::string(continuation_extension), "this is plainly not a state file");
    CHECK_THROWS_WITH(read_continuation_state(foreign), Catch::Matchers::ContainsSubstring("not a continuation state file"));
}

TEST_CASE("ContinuationState::rejects_a_truncated_file") {
    auto rb = make_rigidbody(three_bodies());
    StateFile path;
    write_continuation_state(path, *rb, {});

    std::string contents;
    {
        std::ifstream in(static_cast<const io::File&>(path).str(), std::ios::binary);
        contents.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    REQUIRE(contents.size() > 32);
    {   // keep the header so the magic still matches, then cut the body short
        std::ofstream out(static_cast<const io::File&>(path).str(), std::ios::binary | std::ios::trunc);
        out.write(contents.data(), 32);
    }

    CHECK_THROWS_WITH(read_continuation_state(path), Catch::Matchers::ContainsSubstring("truncated or corrupt"));
}

TEST_CASE("ContinuationState::rejects_a_missing_file") {
    CHECK_THROWS(read_continuation_state(io::File("temp/definitely_not_written.continue")));
}
