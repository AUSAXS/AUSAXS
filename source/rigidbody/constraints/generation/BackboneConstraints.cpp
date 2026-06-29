// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/generation/BackboneConstraints.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

std::vector<std::unique_ptr<IDistanceConstraint>> BackboneConstraints::generate() const {
    if (settings::general::verbose) {console::print_info("\tGenerating backbone constraints for rigid body optimization.");}
    std::vector<std::unique_ptr<IDistanceConstraint>> constraints;

    // Attempt a bond between every pair of bodies and keep the ones that succeed. A DistanceConstraintBond only constructs successfully between backbone-adjacent bodies 
    // (consecutive C-alpha residues within bonding distance), so this naturally selects exactly the backbone bonds regardless of body ordering, and lets 
    // non-adjacent pairs fail silently.
    auto& protein = *manager->molecule;
    for (unsigned int ibody1 = 0; ibody1 < protein.size_body(); ibody1++) {
        for (unsigned int ibody2 = ibody1 + 1; ibody2 < protein.size_body(); ibody2++) {
            try {
                auto constraint = std::make_unique<DistanceConstraintBond>(manager->molecule, ibody1, ibody2);
                if (settings::general::verbose) {
                    std::cout
                        << "\tConstraint created between bodies " << ibody1 << " and " << ibody2 << " on atoms "
                        << form_factor::to_string(constraint->get_atom1().form_factor_type()) << " and " << form_factor::to_string(constraint->get_atom2().form_factor_type())
                    << std::endl;
                }
                constraints.emplace_back(std::move(constraint));
            } catch (const except::base&) {
                // bodies ibody1 and ibody2 are not backbone-adjacent (or too far apart) - skip them
            }
        }
    }

    if (constraints.empty()) {
        console::print_warning("BackboneConstraints::generate: No backbone C-alpha bonds could be identified between any pair of bodies.");
    }

    return constraints;
}
