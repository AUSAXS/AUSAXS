// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/generation/BackboneConstraints.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <utility/Console.h>
#include <utility/Logging.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <string>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

std::vector<std::unique_ptr<IDistanceConstraint>> BackboneConstraints::generate() const {
    console::print_info("\tGenerating backbone constraints for rigid body optimization.");
    std::vector<std::unique_ptr<IDistanceConstraint>> constraints;

    // check every pair of bodies and keep the ones that bond.
    auto& protein = *manager->molecule;
    for (unsigned int ibody1 = 0; ibody1 < protein.size_body(); ibody1++) {
        for (unsigned int ibody2 = ibody1 + 1; ibody2 < protein.size_body(); ibody2++) {
            if (!DistanceConstraintBond::can_bond(manager->molecule, ibody1, ibody2)) {continue;}
            auto constraint = std::make_unique<DistanceConstraintBond>(manager->molecule, ibody1, ibody2);
            logging::log(
                "Constraint created between bodies " + std::to_string(ibody1) + " and " + std::to_string(ibody2) + " on atoms "
                + form_factor::to_string(constraint->get_atom1().form_factor_type()) + " and " + form_factor::to_string(constraint->get_atom2().form_factor_type())
            );
            constraints.emplace_back(std::move(constraint));
        }
    }

    if (constraints.empty()) {
        console::print_warning("BackboneConstraints::generate: No backbone C-alpha bonds could be identified between any pair of bodies.");
        return constraints;
    }
    console::print_text("\tGenerated " + std::to_string(constraints.size()) + " backbone constraints.");

    // A split backbone is a single chain, so a body bonds only to its sequence neighbours: at most two, and never in a loop. The bonds
    // therefore form a forest, for which the number of connected groups is exactly (bodies - bonds) — every bond short of the N-1 a single
    // chain needs leaves one more group behind. A detached group drifts freely, so its position against the data is meaningless; this is
    // reported rather than thrown, since a genuinely multi-chain system is legitimate.
    int nbodies = static_cast<int>(protein.size_body());
    std::vector<unsigned int> links(nbodies, 0);
    for (const auto& constraint : constraints) {
        ++links[constraint->ibody1];
        ++links[constraint->ibody2];
    }

    // more than two bonds means the forest assumption is broken, and with it the group count below
    bool overbonded = false;
    for (int ibody = 0; ibody < nbodies; ++ibody) {
        if (links[ibody] <= 2) {continue;}
        overbonded = true;
        console::print_warning(
            "\tWarning: Body " + std::to_string(ibody) + " is bonded to " + std::to_string(links[ibody]) + " other bodies, but a backbone "
            "chain allows at most two. Check for duplicate or overlapping residue numbering."
        );
    }

    if (int ngroups = nbodies - static_cast<int>(constraints.size()); !overbonded && 1 < ngroups) {
        console::print_warning(
            "\tWarning: The bodies form " + std::to_string(ngroups) + " disconnected groups rather than a single connected structure. Each "
            "group will move independently of the others. The bonds that were formed are listed in the log."
        );
    }

    return constraints;
}
