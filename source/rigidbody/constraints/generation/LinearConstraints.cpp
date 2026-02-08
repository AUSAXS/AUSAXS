// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <settings/GeneralSettings.h>
#include <constants/Constants.h>
#include <utility/Console.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <limits>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

std::vector<std::unique_ptr<IDistanceConstraint>> LinearConstraints::generate() const {
    if (settings::general::verbose) {console::print_info("\tGenerating simple constraints for rigid body optimization.");}
    std::vector<std::unique_ptr<IDistanceConstraint>> constraints;

    auto& protein = *manager->molecule;
    for (unsigned int ibody1 = 0; ibody1 < protein.size_body()-1; ibody1++) {
        unsigned int ibody2 = ibody1 + 1;

        const Body& body1 = protein.get_body(ibody1);
        const Body& body2 = protein.get_body(ibody2);

        double min_dist = std::numeric_limits<double>::max();
        unsigned int min_atom1 = -1, min_atom2 = -1;
        for (unsigned int iatom1 = 0; iatom1 < body1.size_atom(); iatom1++) {
            const AtomFF& atom1 = body1.get_atom(iatom1);
            if (form_factor::to_atom_type(atom1.form_factor_type()) != constants::atom_t::C) {
                continue;
            }

            for (unsigned int iatom2 = 0; iatom2 < body2.size_atom(); iatom2++) {
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

        if (min_atom1 == -1u || min_atom2 == -1u) {
            throw except::unexpected("LinearConstraints::generate: No suitable atoms were found for constraint generation. ");
        }
        constraints.emplace_back(std::make_unique<DistanceConstraintAtom>(manager->molecule, ibody1, ibody2, min_atom1, min_atom2));
        if (settings::general::verbose) {
            std::cout 
                << "\tConstraint created between bodies " << ibody1 << " and " << ibody2 << " on atoms " 
                << form_factor::to_string(body1.get_atom(min_atom1).form_factor_type()) << " and " << form_factor::to_string(body2.get_atom(min_atom2).form_factor_type()) 
            << std::endl;
        }
    }

    if (constraints.empty()) {
        throw except::unexpected("LinearConstraints::generate: No constraints were generated. This is probably a bug.");
    }

    return constraints;
}