// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <settings/RigidBodySettings.h>
#include <constants/Constants.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

double DistanceConstraintCM::evaluate_distance() const {
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

DistanceConstraintCM::DistanceConstraintCM(observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2, std::pair<int, int> isym1, std::pair<int, int> isym2)
    : IDistanceConstraint(molecule, ibody1, ibody2, std::move(isym1), std::move(isym2))
{
    const Body& body1 = molecule->get_body(ibody1);
    const Body& body2 = molecule->get_body(ibody2);

    auto find_cm_atom = [] (const Body& body) -> int {
        double min_distance = std::numeric_limits<double>::max();
        int iatom_cm = -1;
        auto cm = body.get_cm();
        for (unsigned int i = 0; i < body.size_atom(); i++) {
            if (form_factor::to_atom_type(body.get_atom(i).form_factor_type()) != constants::atom_t::C) {continue;}

            double distance = cm.distance2(body.get_atom(i).coordinates());
            if (distance < min_distance) {
                min_distance = distance;
                iatom_cm = i;
            }
        }
        return iatom_cm;
    };

    iatom1 = find_cm_atom(body1);
    iatom2 = find_cm_atom(body2);
    if (iatom1 == -1 || iatom2 == -1) {
        throw except::invalid_argument("DistanceConstraintCM::DistanceConstraintCM: Could not find carbon atoms to represent the center of mass of the two bodies!");
    }
    d_target = body1.get_atom(iatom1).coordinates().distance(body2.get_atom(iatom2).coordinates());
}

double DistanceConstraintCM::evaluate() const {
    return transform(d_target - evaluate_distance());
}

double DistanceConstraintCM::transform(double offset) {
    return offset*offset*offset*offset*10;
}