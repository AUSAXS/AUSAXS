// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintFunctions.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <array>
#include <limits>
#include <optional>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;

namespace {
    // Indices of the first and last C-alpha atoms of a body, i.e. its terminal backbone atoms. A backbone bond can only ever attach at these, 
    // so they are the only selection candidates. Returns {-1, -1} if the body contains no C-alpha atoms.
    std::pair<int, int> terminal_calphas(const data::Body& body) {
        const auto& backbone = *body.get_metadata()->backbone;
        int first = -1, last = -1;
        for (int i = 0; i < static_cast<int>(body.size_atom()); ++i) {
            if (backbone[i] != data::backbone_t::c_alpha) {continue;}
            if (first == -1) {first = i;}
            last = i;
        }
        return {first, last};
    }

    // Residue sequence id of atom i, if sequence metadata is available for this body.
    std::optional<int> residue_seq(const data::Body& body, unsigned int i) {
        const auto& md = body.get_metadata();
        if (md && md->residue_seq) {return (*md->residue_seq)[i];}
        return std::nullopt;
    }
}

DistanceConstraintBond::DistanceConstraintBond(observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2)
    : DistanceConstraintAtom(molecule, ibody1, ibody2) {
    const data::Body& body1 = molecule->get_body(ibody1);
    const data::Body& body2 = molecule->get_body(ibody2);

    assert((body1.get_metadata() && body1.get_metadata()->backbone) && "Backbone metadata must be present for DistanceConstraintBond to work.");
    assert((body2.get_metadata() && body2.get_metadata()->backbone) && "Backbone metadata must be present for DistanceConstraintBond to work.");

    // Only the terminal C-alphas of each body can form a backbone bond, so the candidates are the
    // (up to) four combinations of their endpoints.
    auto [first1, last1] = terminal_calphas(body1);
    auto [first2, last2] = terminal_calphas(body2);
    if (first1 == -1 || first2 == -1) {
        throw except::invalid_argument("DistanceConstraintBond::DistanceConstraintBond: Could not find suitable atoms to represent the bond between the two bodies!");
    }
    const std::array<std::pair<int, int>, 4> candidates = {{ {first1, first2}, {first1, last2}, {last1, first2}, {last1, last2} }};

    // When residue metadata is available we require a *sequential* pair (consecutive residue ids): only a real peptide bond is a backbone bond. 
    // Without that metadata we fall back to the closest terminal pair, so the constraint still works for synthetic/in-memory molecules.
    bool have_seq = residue_seq(body1, first1).has_value() && residue_seq(body2, first2).has_value();
    double best_distance = std::numeric_limits<double>::max();
    iatom1 = -1; iatom2 = -1;
    for (auto [i, j] : candidates) {
        double distance = evaluate_distance(i, j);
        bool sequential = false;
        if (have_seq) {
            int s1 = *residue_seq(body1, i), s2 = *residue_seq(body2, j);
            sequential = (s1 - s2 == 1 || s2 - s1 == 1);
        }
        if ((have_seq && !sequential) || distance >= best_distance) {continue;}
        best_distance = distance;
        iatom1 = i;
        iatom2 = j;
    }
    if (iatom1 == -1 || iatom2 == -1) {
        throw except::invalid_argument("DistanceConstraintBond::DistanceConstraintBond: The two bodies are not backbone-adjacent; no sequential C-alpha pair could be found!");
    }

    d_target = best_distance;
    if (d_target > 4) {
        auto& atom1 = body1.get_atom(iatom1);
        auto& atom2 = body2.get_atom(iatom2);
        throw except::invalid_argument(
            "DistanceConstraint::DistanceConstraint: The atoms being constrained are too far apart!\n"
            "Atom 1: " + form_factor::to_string(atom1.form_factor_type()) + " in body " + std::to_string(ibody1) + " at " + atom1.coordinates().to_string() + "\n"
            "Atom 2: " + form_factor::to_string(atom2.form_factor_type()) + " in body " + std::to_string(ibody2) + " at " + atom2.coordinates().to_string() + "\n"
        );
    }
}

double DistanceConstraintBond::evaluate() const {
    double distance = evaluate_distance(iatom1, iatom2);
    double offset = distance - d_target;
    return functions::between_atoms(offset);
}