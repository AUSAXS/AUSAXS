// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/DistanceConstraintFunctions.h>
#include <constants/Constants.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

double DistanceConstraintAtom::evaluate_distance() const {
    auto atom1 = molecule->get_body(ibody1).get_atom(iatom1).coordinates();
    auto atom2 = molecule->get_body(ibody2).get_atom(iatom2).coordinates();
    if (isym1.first != -1) {
        const auto& sym = molecule->get_body(ibody1).symmetry().get(isym1.first);
        auto transform = sym.get_transform<double>(atom1, isym1.second);
        atom1 = transform(atom1);
    }
    if (isym2.first != -1) {
        const auto& sym = molecule->get_body(ibody2).symmetry().get(isym2.first);
        auto transform = sym.get_transform<double>(atom2, isym2.second);
        atom2 = transform(atom2);
    }
    return atom1.distance(atom2);
}

DistanceConstraintAtom::DistanceConstraintAtom(
    observer_ptr<const data::Molecule> molecule, int ibody1, int iatom1, int ibody2, int iatom2, std::pair<int, int> isym1, std::pair<int, int> isym2
) : IDistanceConstraint(molecule, ibody1, ibody2, std::move(isym1), std::move(isym2)) {
    this->iatom1 = iatom1;
    this->iatom2 = iatom2;
    
    if (iatom1 < 0 || static_cast<size_t>(iatom1) >= molecule->get_body(ibody1).size_atom()) {
        throw except::invalid_argument("DistanceConstraintAtom::DistanceConstraintAtom: Invalid atom index " + std::to_string(iatom1) + " for body " + std::to_string(ibody1));
    }
    if (iatom2 < 0 || static_cast<size_t>(iatom2) >= molecule->get_body(ibody2).size_atom()) {
        throw except::invalid_argument("DistanceConstraintAtom::DistanceConstraintAtom: Invalid atom index " + std::to_string(iatom2) + " for body " + std::to_string(ibody2));
    }

    const auto& atom1 = molecule->get_body(ibody1).get_atom(iatom1);
    const auto& atom2 = molecule->get_body(ibody2).get_atom(iatom2);

    // Only allow constraints between carbon atoms (backbone C-alpha)
    if (form_factor::to_atom_type(atom1.form_factor_type()) != constants::atom_t::C ||
        form_factor::to_atom_type(atom2.form_factor_type()) != constants::atom_t::C) {
        throw except::invalid_argument("DistanceConstraintAtom::DistanceConstraintAtom: Constraints only make sense between carbon atoms of the backbone!");
    }

    // Constraints within the same body don't make sense
    if (molecule->get_body(ibody1).get_uid() == molecule->get_body(ibody2).get_uid()) {
        throw except::invalid_argument("DistanceConstraintAtom::DistanceConstraintAtom: Cannot create a constraint between atoms in the same body!");
    }

    d_target = evaluate_distance();
}

DistanceConstraintAtom::DistanceConstraintAtom(
    observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2,
    std::pair<int, int> isym1, std::pair<int, int> isym2
) : IDistanceConstraint(molecule, ibody1, ibody2, std::move(isym1), std::move(isym2)) {}

DistanceConstraintAtom::DistanceConstraintAtom(
    observer_ptr<const data::Molecule> molecule,
    const data::AtomFF& atom1,
    const data::AtomFF& atom2
) : IDistanceConstraint(molecule, -1, -1, {-1, -1}, {-1, -1}) {
    // Only allow constraints between carbon atoms (backbone C-alpha)
    if (form_factor::to_atom_type(atom1.form_factor_type()) != constants::atom_t::C ||
        form_factor::to_atom_type(atom2.form_factor_type()) != constants::atom_t::C) {
        throw except::invalid_argument("DistanceConstraintAtom::DistanceConstraintAtom: Constraints only make sense between carbon atoms of the backbone!");
    }

    // Find which bodies contain these atoms
    int found_ibody1 = -1, found_ibody2 = -1;
    int found_iatom1 = -1, found_iatom2 = -1;
    for (unsigned int ibody = 0; ibody < molecule->size_body(); ++ibody) {
        const auto& body = molecule->get_body(ibody);
        for (unsigned int iatom = 0; iatom < body.size_atom(); ++iatom) {
            const auto& atom = body.get_atom(iatom);
            if (atom1 == atom) {
                found_ibody1 = ibody;
                found_iatom1 = iatom;
                break; // atoms must be from different bodies
            } else if (atom2 == atom) {
                found_ibody2 = ibody;
                found_iatom2 = iatom;
                break;
            }
        }
    }

    if (found_ibody1 == -1 || found_ibody2 == -1) {
        throw except::invalid_argument("DistanceConstraintAtom::DistanceConstraintAtom: Could not determine host bodies for the two atoms.");
    }

    // Constraints within the same body don't make sense
    if (found_ibody1 == found_ibody2) {
        throw except::invalid_argument("DistanceConstraintAtom::DistanceConstraintAtom: Cannot create a constraint between atoms in the same body!");
    }

    ibody1 = found_ibody1;
    ibody2 = found_ibody2;
    iatom1 = found_iatom1;
    iatom2 = found_iatom2;

    d_target = atom1.coordinates().distance(atom2.coordinates());
}

double DistanceConstraintAtom::evaluate() const {
    assert(d_target != 0);
    return functions::between_atoms(d_target - evaluate_distance());
}