// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/generation/LinearConstraints.h>
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

std::vector<std::unique_ptr<IDistanceConstraint>> LinearConstraints::generate() const {
    if (settings::general::verbose) {console::print_info("\tGenerating simple constraints for rigid body optimization.");}
    std::vector<std::unique_ptr<IDistanceConstraint>> constraints;

    // A linear chain simply places one bond constraint between each pair of sequential bodies.
    // DistanceConstraintBond already picks the most appropriate C-alpha pair between two bodies,
    // so we delegate the per-pair atom selection to it instead of duplicating that logic here.
    auto& protein = *manager->molecule;
    for (unsigned int ibody1 = 0; ibody1 + 1 < protein.size_body(); ibody1++) {
        unsigned int ibody2 = ibody1 + 1;
        auto constraint = std::make_unique<DistanceConstraintBond>(manager->molecule, ibody1, ibody2);
        if (settings::general::verbose) {
            std::cout
                << "\tConstraint created between bodies " << ibody1 << " and " << ibody2 << " on atoms "
                << form_factor::to_string(constraint->get_atom1().form_factor_type()) << " and " << form_factor::to_string(constraint->get_atom2().form_factor_type())
            << std::endl;
        }
        constraints.emplace_back(std::move(constraint));
    }

    if (constraints.empty()) {
        throw except::unexpected("LinearConstraints::generate: No constraints were generated. This is probably a bug.");
    }

    return constraints;
}