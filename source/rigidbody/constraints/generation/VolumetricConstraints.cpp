// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/generation/VolumetricConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <settings/RigidBodySettings.h>
#include <settings/GeneralSettings.h>
#include <constants/Constants.h>
#include <utility/Console.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <limits>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

std::vector<DistanceConstraint> VolumetricConstraints::generate() const {
    if (settings::general::verbose) {console::print_info("\tGenerating simple constraints for rigid body optimization.");}
    std::vector<DistanceConstraint> constraints;

    auto& protein = *manager->molecule;
    for (unsigned int ibody1 = 0; ibody1 < protein.get_bodies().size(); ibody1++) {
        for (unsigned int ibody2 = ibody1+1; ibody2 < protein.get_bodies().size(); ibody2++) {
            const Body& body1 = protein.get_body(ibody1);
            const Body& body2 = protein.get_body(ibody2);

            double min_dist = std::numeric_limits<double>::max();
            int min_atom1 = -1, min_atom2 = -1;
            for (unsigned int iatom1 = 0; iatom1 < body1.get_atoms().size(); iatom1++) {
                const AtomFF& atom1 = body1.get_atom(iatom1);
                if (form_factor::to_atom_type(atom1.form_factor_type()) != constants::atom_t::C) {
                    continue;
                }

                for (unsigned int iatom2 = 0; iatom2 < body2.get_atoms().size(); iatom2++) {
                    const AtomFF& atom2 = body2.get_atom(iatom2);
                    if (form_factor::to_atom_type(atom2.form_factor_type()) != constants::atom_t::C) {
                        continue;
                    }

                    double dist = atom1.coordinates().distance(atom2.coordinates());
                    if (dist > min_dist) {continue;}

                    min_dist = dist;
                    min_atom1 = iatom1;
                    min_atom2 = iatom2;
                }
            }

            // no carbon atoms found
            if (min_atom1 == -1 || min_atom2 == -1) {continue;}

            // check if the bodies are close enough for a constraint to make sense
            if (min_dist > settings::rigidbody::bond_distance) {continue;} 

            constraints.emplace_back(manager->molecule, ibody1, ibody2, min_atom1, min_atom2);
            if (settings::general::verbose) {
                std::cout 
                    << "\tConstraint created between bodies " << ibody1 << " and " << ibody2 << " on atoms " 
                    << form_factor::to_string(body1.get_atom(min_atom1).form_factor_type()) << " and " << form_factor::to_string(body2.get_atom(min_atom2).form_factor_type()) 
                << std::endl;
            }
        }
    }

    if (constraints.empty()) {
        throw except::unexpected("rigidbody::constraints::generation::VolumeConstraints::generate: No constraints were generated. This is probably a bug.");
    }

    return constraints;
}