// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintFunctions.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;

namespace {
    bool is_constraint_candidate(const data::Body& body, unsigned int i) {
        const auto& md = body.get_metadata();
        return (*md->backbone)[i] == data::backbone_t::c_alpha;
    }

    // Residue sequence id of atom i, if sequence metadata is available for this body.
    std::optional<int> residue_seq(const data::Body& body, unsigned int i) {
        const auto& md = body.get_metadata();
        if (md && md->residue_seq) {return (*md->residue_seq)[i];}
        return std::nullopt;
    }
}

DistanceConstraintBond::DistanceConstraintBond(observer_ptr<const data::Molecule> molecule,
    int ibody1, int ibody2, std::pair<int, int> isym1, std::pair<int, int> isym2
) : DistanceConstraintAtom(molecule, ibody1, ibody2, std::move(isym1), std::move(isym2)) {
    const data::Body& body1 = molecule->get_body(ibody1);
    const data::Body& body2 = molecule->get_body(ibody2);

    assert((body1.get_metadata() && body1.get_metadata()->backbone) && "Backbone metadata must be present for DistanceConstraintBond to work.");
    assert((body2.get_metadata() && body2.get_metadata()->backbone) && "Backbone metadata must be present for DistanceConstraintBond to work.");

    double min_distance = std::numeric_limits<double>::max();     // closest pair overall
    double min_seq_distance = std::numeric_limits<double>::max(); // closest pair of sequential residues
    int seq_atom1 = -1, seq_atom2 = -1;
    for (int i = 0; i < static_cast<int>(body1.size_atom()); ++i) {
        if (!is_constraint_candidate(body1, i)) {
            continue;
        }
        std::optional<int> seq1 = residue_seq(body1, i);

        for (int j = 0; j < static_cast<int>(body2.size_atom()); ++j) {
            if (!is_constraint_candidate(body2, j)) {
                continue;
            }

            double distance = evaluate_distance(i, j);
            if (distance < min_distance) {
                min_distance = distance;
                iatom1 = i;
                iatom2 = j;
            }

            std::optional<int> seq2 = residue_seq(body2, j);
            if (seq1 && seq2 && (*seq1 - *seq2 == 1 || *seq2 - *seq1 == 1) && distance < min_seq_distance) {
                min_seq_distance = distance;
                seq_atom1 = i;
                seq_atom2 = j;
            }
        }
    }
    if (seq_atom1 != -1) { // a sequential pair was found; prefer it over the merely-closest pair
        iatom1 = seq_atom1;
        iatom2 = seq_atom2;
        min_distance = min_seq_distance;
    }
    if (iatom1 == -1 || iatom2 == -1) {
        throw except::invalid_argument("DistanceConstraintBond::DistanceConstraintBond: Could not find suitable atoms to represent the bond between the two bodies!");
    }

    d_target = min_distance;
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