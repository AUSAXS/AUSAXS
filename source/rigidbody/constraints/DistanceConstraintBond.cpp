// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintFunctions.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs::rigidbody::constraints;

DistanceConstraintBond::DistanceConstraintBond(observer_ptr<const data::Molecule> molecule,
    int ibody1, int ibody2, std::pair<int, int> isym1, std::pair<int, int> isym2
) : DistanceConstraintAtom(molecule, ibody1, ibody2, std::move(isym1), std::move(isym2)) {
    const data::Body& body1 = molecule->get_body(ibody1);
    const data::Body& body2 = molecule->get_body(ibody2);

    // find the closest atoms in the two bodies
    double min_distance = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < body1.size_atom(); i++) {
        if (form_factor::to_atom_type(body1.get_atom(i).form_factor_type()) != constants::atom_t::C) {
            continue;
        }

        for (unsigned int j = 0; j < body2.size_atom(); j++) {
            if (form_factor::to_atom_type(body2.get_atom(j).form_factor_type()) != constants::atom_t::C) {
                continue;
            }

            double distance = body1.get_atom(i).coordinates().distance2(body2.get_atom(j).coordinates());
            if (distance < min_distance) {
                min_distance = distance;
                iatom1 = i;
                iatom2 = j;
            }
        }
    }
    if (iatom1 == -1 || iatom2 == -1) {
        throw except::invalid_argument("DistanceConstraintBond::DistanceConstraintBond: Could not find carbon atoms to represent the bond between the two bodies!");
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
    double distance = evaluate_distance();
    if (d_target < distance) {return 0;}
    double offset = distance - d_target;
    return functions::attractor_repulsor(offset);
}