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
    d_target = evaluate_distance();
}

DistanceConstraintAtom::DistanceConstraintAtom(
    observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2,
    std::pair<int, int> isym1, std::pair<int, int> isym2
) : IDistanceConstraint(molecule, ibody1, ibody2, std::move(isym1), std::move(isym2)) {}

double DistanceConstraintAtom::evaluate() const {
    assert(d_target != 0);
    return functions::between_atoms(d_target - evaluate_distance());
}