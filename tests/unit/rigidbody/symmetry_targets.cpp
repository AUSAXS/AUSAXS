// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/selection/SymmetryTargets.h>
#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/selection/ManualSelect.h>
#include <rigidbody/selection/ParameterMask.h>
#include <rigidbody/selection/ParameterMaskStrategy.h>
#include <rigidbody/Rigidbody.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

#include <memory>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::selection;

namespace {
    // Four single-atom bodies, no constraints. Symmetries are attached per test.
    std::unique_ptr<Rigidbody> make_rigidbody(int nbodies = 4) {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;

        std::vector<Body> bodies;
        for (int i = 0; i < nbodies; ++i) {
            bodies.emplace_back(std::vector{AtomFF({5.0*i, 0, 0}, form_factor::form_factor_t::C)});
        }
        return std::make_unique<Rigidbody>(Molecule{std::move(bodies)});
    }

    // Tie every body into one shared symmetry, exactly as SplitElement does: body 0 owns it, the rest hold non-owning views.
    void share_symmetry(Rigidbody& rb, symmetry::type base) {
        std::vector<int> participants(rb.molecule.size_body());
        for (int i = 0; i < static_cast<int>(participants.size()); ++i) {participants[i] = i;}

        rb.molecule.get_body(0).symmetry().add(std::make_unique<symmetry::ReferenceSymmetry>(
            symmetry::get(base), participants, std::vector<int>(participants.size(), 0), &rb.molecule
        ));
        for (std::size_t b = 1; b < rb.molecule.size_body(); ++b) {
            rb.molecule.get_body(b).symmetry().add(std::make_unique<symmetry::ReferenceSymmetryView>(&rb.molecule, 0, 0));
        }
        rb.symmetry_targets->invalidate();
    }
}

TEST_CASE("symmetry::is_optimizable sees through the reference wrapper") {
    auto rb = make_rigidbody(2);

    SECTION("a plain symmetry is drivable") {
        CHECK(symmetry::is_optimizable(*symmetry::get(symmetry::type::c2)));
        CHECK(symmetry::is_optimizable(*symmetry::create("p2-c3")));
    }

    SECTION("a ReferenceSymmetry is judged by the base it wraps, not by being a wrapper") {
        symmetry::ReferenceSymmetry ref(symmetry::get(symmetry::type::c3), {0, 1}, {0, 0}, &rb->molecule);
        CHECK(symmetry::is_optimizable(ref));
    }

    SECTION("a ReferenceSymmetryView is not drivable: it exposes no parameters of its own") {
        symmetry::ReferenceSymmetryView view(&rb->molecule, 0, 0);
        CHECK_FALSE(symmetry::is_optimizable(view));
    }
}

TEST_CASE("SymmetryTargets: only the owning slot of a shared symmetry is drivable") {
    auto rb = make_rigidbody(4);
    REQUIRE(rb->symmetry_targets->empty());

    share_symmetry(*rb, symmetry::type::c2);

    // every body declares one symmetry, but only body 0 owns it
    for (std::size_t b = 0; b < rb->molecule.size_body(); ++b) {
        REQUIRE(rb->molecule.get_body(b).size_symmetry() == 1);
    }
    REQUIRE(rb->symmetry_targets->size() == 1);
    CHECK(rb->symmetry_targets->all()[0].ibody == 0);
    CHECK(rb->symmetry_targets->all()[0].isymmetry == 0);

    CHECK(rb->symmetry_targets->body_targets(0) == std::vector<unsigned int>{0});
    for (unsigned int b = 1; b < rb->molecule.size_body(); ++b) {
        CHECK(rb->symmetry_targets->body_targets(b).empty());
    }
}

TEST_CASE("SymmetryTargets: the pool is rebuilt after invalidation") {
    auto rb = make_rigidbody(2);
    REQUIRE(rb->symmetry_targets->empty());

    // reading the pool caches it; a subsequent change must not be visible until it is invalidated
    rb->molecule.get_body(0).symmetry().add(symmetry::type::c2);
    CHECK(rb->symmetry_targets->empty());

    rb->symmetry_targets->invalidate();
    REQUIRE(rb->symmetry_targets->size() == 1);
    CHECK(rb->symmetry_targets->all()[0].ibody == 0);

    rb->molecule.get_body(1).symmetry().add(symmetry::type::c3);
    rb->symmetry_targets->invalidate();
    REQUIRE(rb->symmetry_targets->size() == 2);
    CHECK(rb->symmetry_targets->all()[1].ibody == 1);
}

TEST_CASE("SelectionStrategies: a symmetry-only mask never selects a body that cannot move") {
    auto rb = make_rigidbody(4);
    share_symmetry(*rb, symmetry::type::c2);
    REQUIRE(rb->symmetry_targets->size() == 1); // bodies 1-3 hold views only

    auto check = [&](BodySelectStrategy& selector) {
        for (int i = 0; i < 50; ++i) {
            auto target = selector.next(ParameterMask::symmetry_only());
            CHECK(target.ibody == 0);                 // the only body owning a drivable slot
            CHECK(target.isymmetry == 0);
            CHECK(target.iconstraint == -1);          // driving a symmetry moves no atoms of the host body
        }
    };

    SECTION("RandomBodySelect")     {RandomBodySelect s(rb.get());     check(s);}
    SECTION("SequentialBodySelect") {SequentialBodySelect s(rb.get()); check(s);}

    SECTION("a real-space mask still draws freely over all bodies") {
        RandomBodySelect selector(rb.get());
        bool saw_view_body = false;
        for (int i = 0; i < 50; ++i) {
            auto target = selector.next(ParameterMask::real_only());
            CHECK(target.isymmetry == -1);
            saw_view_body |= (target.ibody != 0);
        }
        INFO("bodies holding only shared-symmetry views are still valid pose targets");
        CHECK(saw_view_body);
    }
}

TEST_CASE("SelectionStrategies: next_mask carries the chosen slot into the mask") {
    auto rb = make_rigidbody(2);
    rb->molecule.get_body(1).symmetry().add(symmetry::type::c2);
    rb->molecule.get_body(1).symmetry().add(symmetry::type::c3);
    rb->symmetry_targets->invalidate();

    ManualSelect selector(rb.get(), 1, 1); // the second declared symmetry of body 1
    selector.set_mask_strategy(std::make_unique<SymmetryOnlyMaskStrategy>());

    auto selection = selector.next_mask();
    CHECK(selection.ibody == 1);
    REQUIRE(selection.mask.target_symmetry.has_value());
    CHECK(selection.mask.target_symmetry.value() == 1);
}

TEST_CASE("SymmetryTargets::resolve maps a shared symmetry onto the body owning it") {
    auto rb = make_rigidbody(3);
    share_symmetry(*rb, symmetry::type::c2);

    SECTION("the owner resolves to itself") {
        auto target = rb->symmetry_targets->resolve(0, 0);
        REQUIRE(target.has_value());
        CHECK(target->ibody == 0);
        CHECK(target->isymmetry == 0);
    }

    SECTION("every participant resolves to the owner, so any of their names reaches the same parameters") {
        for (unsigned int b = 1; b < rb->molecule.size_body(); ++b) {
            auto target = rb->symmetry_targets->resolve(b, 0);
            REQUIRE(target.has_value());
            CHECK(target->ibody == 0);
            CHECK(target->isymmetry == 0);
        }
    }

    SECTION("a slot that does not exist resolves to nothing") {
        CHECK_FALSE(rb->symmetry_targets->resolve(0, 5).has_value());
        CHECK_FALSE(rb->symmetry_targets->resolve(99, 0).has_value());
    }
}

TEST_CASE("ManualSelect: naming a shared symmetry through any participant drives the owner's copy") {
    auto rb = make_rigidbody(3);
    share_symmetry(*rb, symmetry::type::c2);

    // b2s1 / b2s1r2 resolve to {body 1, symmetry 0}, which is only a view; the parameters live on body 0
    auto named_through = GENERATE(0u, 1u, 2u);
    ManualSelect selector(rb.get(), named_through, 0);
    selector.set_mask_strategy(std::make_unique<SymmetryOnlyMaskStrategy>());

    auto selection = selector.next_mask();
    CHECK(selection.ibody == 0);
    CHECK(selection.iconstraint == -1);
    REQUIRE(selection.mask.target_symmetry.has_value());
    CHECK(selection.mask.target_symmetry.value() == 0);
}

TEST_CASE("SelectionStrategies: an undrivable target is rejected rather than silently doing nothing") {
    SECTION("a slot backed by nothing at all throws") {
        auto rb = make_rigidbody(2);
        rb->molecule.get_body(0).symmetry().add(symmetry::type::c2);
        rb->symmetry_targets->invalidate();
        CHECK_NOTHROW(ManualSelect(rb.get(), 0, 0));
        CHECK_THROWS(ManualSelect(rb.get(), 1, 0)); // body 1 declares no symmetry at all
    }

    SECTION("a symmetry-only run over a molecule with no drivable symmetry throws") {
        auto bare = make_rigidbody(3);
        RandomBodySelect selector(bare.get());
        selector.set_mask_strategy(std::make_unique<SymmetryOnlyMaskStrategy>());
        CHECK_THROWS(selector.next_mask());
    }
}

TEST_CASE("ManualSelect: selecting a body honours the mask it is given") {
    // the molecule as a whole has a drivable slot - body 0 owns the shared symmetry - so next_mask's molecule-wide guard stays silent and the selected body
    // is the only thing standing between the optimizer and a run that perturbs nothing
    auto rb = make_rigidbody(3);
    share_symmetry(*rb, symmetry::type::c2);
    REQUIRE(rb->symmetry_targets->size() == 1);

    SECTION("a body holding only a shared-symmetry view cannot satisfy a symmetry-only mask") {
        ManualSelect selector(rb.get(), 1); // declares one symmetry, but only as a view onto body 0's copy
        selector.set_mask_strategy(std::make_unique<SymmetryOnlyMaskStrategy>());
        CHECK_THROWS(selector.next_mask());
    }

    SECTION("neither can a body declaring no symmetry at all") {
        auto bare = make_rigidbody(2);
        bare->molecule.get_body(0).symmetry().add(symmetry::type::c2);
        bare->symmetry_targets->invalidate();

        ManualSelect selector(bare.get(), 1);
        selector.set_mask_strategy(std::make_unique<SymmetryOnlyMaskStrategy>());
        CHECK_THROWS(selector.next_mask());
    }

    SECTION("a body that does own a drivable slot moves all of its symmetries together") {
        ManualSelect selector(rb.get(), 0);
        selector.set_mask_strategy(std::make_unique<SymmetryOnlyMaskStrategy>());

        auto selection = selector.next_mask();
        CHECK(selection.ibody == 0);
        // selecting the body rather than one of its slots leaves the target unset, so every symmetry it declares moves alike
        CHECK_FALSE(selection.mask.target_symmetry.has_value());
    }

    SECTION("a mask that leaves the pose free is satisfied by any body") {
        ManualSelect selector(rb.get(), 1);
        selector.set_mask_strategy(std::make_unique<AllMaskStrategy>());
        CHECK_NOTHROW(selector.next_mask());
    }
}
