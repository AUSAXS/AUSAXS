// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/atoms/AtomMetadata.h>
#include <data/Molecule.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <settings/All.h>

#include <support/rb_metadata.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

struct fixture {
    AtomFF a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3 = AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
    AtomFF a4 = AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
    AtomFF a5 = AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a8 = AtomFF({ 1,  1,  1}, form_factor::form_factor_t::NH);

    Body b1 = Body(std::vector<AtomFF>{a1, a2});
    Body b2 = Body(std::vector<AtomFF>{a3, a4});
    Body b3 = Body(std::vector<AtomFF>{a5, a8}); // b3 has only one C atom
    std::vector<Body> ap = {b1, b2, b3};

    fixture() {test::mark_backbone_carbons(ap);}
};

TEST_CASE_METHOD(fixture, "DistanceConstraintBond::constructor") {
    settings::molecule::implicit_hydrogens = false;
    Molecule protein = Molecule(ap);

    SECTION("selects the closest C-C pair between the two bodies") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);
        // a1(-1,-1,-1) to a3(-1,-1,1): dist=2 is the closest pair (found first)
        CHECK(c.ibody1 == 0);
        CHECK(c.ibody2 == 1);
        CHECK(c.get_atom1() == a1);
        CHECK(c.get_atom2() == a3);
        CHECK_THAT(c.d_target, Catch::Matchers::WithinAbs(2.0, 1e-9));
    }

    SECTION("throws if atoms are too far apart") {
        Body far1 = Body(std::vector<AtomFF>{AtomFF({0,  0, 0}, form_factor::form_factor_t::C)});
        Body far2 = Body(std::vector<AtomFF>{AtomFF({0, 10, 0}, form_factor::form_factor_t::C)});
        std::vector<Body> far_ap = {far1, far2};
        Molecule far_protein(far_ap);
        test::mark_backbone_carbons(far_protein);
        CHECK_THROWS(constraints::DistanceConstraintBond(&far_protein, 0, 1));
    }
}

TEST_CASE_METHOD(fixture, "DistanceConstraintBond::evaluate") {
    settings::molecule::implicit_hydrogens = false;
    Molecule protein = Molecule(ap);

    SECTION("relaxed") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);
        CHECK(c.evaluate() == 0);
    }

    SECTION("stretched - penalty above d_target") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);

        // Move body0 away from body1 along z: distance increases beyond d_target
        protein.get_body(0).translate(Vector3<double>(0, 0, -2));
        CHECK(c.evaluate() != 0);
    }

    SECTION("compressed - penalty below d_target") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);

        // Move body0 towards body1 along z: distance decreases below d_target
        protein.get_body(0).translate(Vector3<double>(0, 0, 1));
        CHECK(c.evaluate() != 0);

        protein.get_body(0).translate(Vector3<double>(0, 0, -1));
        CHECK(c.evaluate() == 0);
    }

    SECTION("unrelated body translation does not affect constraint") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);
        protein.get_body(2).translate(Vector3<double>(5, 0, 0));
        CHECK(c.evaluate() == 0);
    }

    SECTION("symmetry on unrelated body does not affect constraint") {
        constraints::DistanceConstraintBond c(&protein, 0, 1);
        protein.get_body(2).symmetry().add(symmetry::type::c2);
        auto* sym2 = static_cast<symmetry::CyclicSymmetry*>(protein.get_body(2).symmetry().get(0));
        sym2->_initial_relation.translation = {0.5, 0, 0};
        CHECK(c.evaluate() == 0);
    }
}

TEST_CASE("DistanceConstraintBond: prefers sequential C-alpha pairs") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    // This structure used to produce bond constraints between non-sequential C-alpha atoms (the
    // geometrically closest pair) instead of the sequential pair joining the two bodies. The generated
    // bond constraints identify that pair through the C-alpha and residue-sequence metadata.
    Rigidbody rigidbody = BodySplitter::split("tests/files/LAR1-4.pdb", {9, 99, 202, 292});
    rigidbody.constraints->generate_constraints(settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone);

    const auto& constraints = rigidbody.constraints->discoverable_constraints;
    REQUIRE(constraints.size() == 4); // one bond between each of the 5 sequential bodies

    for (const auto& c : constraints) {
        const auto& md1 = c->get_body1().get_metadata();
        const auto& md2 = c->get_body2().get_metadata();
        REQUIRE((md1 && md1->backbone && md1->residue_seq));
        REQUIRE((md2 && md2->backbone && md2->residue_seq));

        // both constrained atoms must be C-alpha backbone atoms...
        CHECK((*md1->backbone)[c->iatom1] == data::backbone_t::c_alpha);
        CHECK((*md2->backbone)[c->iatom2] == data::backbone_t::c_alpha);

        // ...belonging to sequential residues (adjacent residue sequence ids)
        int seq1 = (*md1->residue_seq)[c->iatom1];
        int seq2 = (*md2->residue_seq)[c->iatom2];
        INFO("constraint between residues " << seq1 << " and " << seq2);
        CHECK((seq1 - seq2 == 1 || seq2 - seq1 == 1));
    }
}

// Build two single-atom bodies carrying full backbone/residue/chain metadata, separated along x.
static std::vector<Body> two_residues(int seq1, char chain1, int seq2, char chain2, double separation) {
    auto make = [] (Vector3<double> pos, int seq, char chain) {
        Body b(std::vector<AtomFF>{AtomFF(pos, form_factor::form_factor_t::C)});
        AtomMetadata md;
        md.backbone.emplace(1, data::backbone_t::c_alpha);
        md.residue_seq.emplace(1, seq);
        md.chain_id.emplace(1, chain);
        b.set_metadata(std::move(md));
        return b;
    };
    return {make({0, 0, 0}, seq1, chain1), make({separation, 0, 0}, seq2, chain2)};
}

TEST_CASE("DistanceConstraintBond: residues missing from the model") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    // A deposited model routinely omits disordered residues. The two fragments on either side of such a gap are still
    // covalently connected, so they must still bond - the gap is paid for in the distance the bond is allowed to span.
    SECTION("one missing residue still bonds") {
        auto ap = two_residues(10, 'A', 12, 'A', 6);
        Molecule protein(ap);
        REQUIRE(constraints::DistanceConstraintBond::can_bond(&protein, 0, 1));
        constraints::DistanceConstraintBond c(&protein, 0, 1);
        CHECK_THAT(c.d_target, Catch::Matchers::WithinAbs(6.0, 1e-9));
    }

    SECTION("four missing residues still bond") {
        auto ap = two_residues(10, 'A', 15, 'A', 16);
        Molecule protein(ap);
        CHECK(constraints::DistanceConstraintBond::can_bond(&protein, 0, 1));
    }

    // ...but only as far as the missing stretch can physically reach: one residue step spans at most ~4A.
    SECTION("a gap the missing residues cannot span does not bond") {
        auto ap = two_residues(10, 'A', 12, 'A', 9);
        Molecule protein(ap);
        CHECK_FALSE(constraints::DistanceConstraintBond::can_bond(&protein, 0, 1));
        CHECK_THROWS(constraints::DistanceConstraintBond(&protein, 0, 1));
    }

    SECTION("a stretch too long to be a single connection does not bond") {
        auto ap = two_residues(10, 'A', 40, 'A', 3.8);
        Molecule protein(ap);
        CHECK_FALSE(constraints::DistanceConstraintBond::can_bond(&protein, 0, 1));
    }

    // Residue ids restart on every chain, so two chains regularly carry consecutive ids without being connected at all.
    SECTION("adjacent residue ids on different chains do not bond") {
        auto ap = two_residues(10, 'A', 11, 'B', 3.8);
        Molecule protein(ap);
        CHECK_FALSE(constraints::DistanceConstraintBond::can_bond(&protein, 0, 1));

        auto same = two_residues(10, 'A', 11, 'A', 3.8);
        Molecule same_protein(same);
        CHECK(constraints::DistanceConstraintBond::can_bond(&same_protein, 0, 1));
    }
}
