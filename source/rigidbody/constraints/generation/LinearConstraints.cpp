#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <data/Protein.h>

using namespace rigidbody;

std::vector<std::shared_ptr<rigidbody::DistanceConstraint>> LinearConstraints::generate() const {
    if (setting::general::verbose) {utility::print_info("\tGenerating simple constraints for rigid body optimization.");}
    std::vector<std::shared_ptr<rigidbody::DistanceConstraint>> constraints;

    auto& protein = *manager->protein;
    for (unsigned int ibody1 = 0; ibody1 < protein.bodies.size()-1; ibody1++) {
        unsigned int ibody2 = ibody1 + 1;

        const Body& body1 = protein.body(ibody1);
        const Body& body2 = protein.body(ibody2);

        double min_dist = std::numeric_limits<double>::max();
        int min_atom1 = -1, min_atom2 = -1;
        for (unsigned int iatom1 = 0; iatom1 < body1.atoms().size(); iatom1++) {
            const Atom& atom1 = body1.atoms(iatom1);
            if (atom1.element != constants::symbols::carbon) {continue;}

            for (unsigned int iatom2 = 0; iatom2 < body2.atoms().size(); iatom2++) {
                const Atom& atom2 = body2.atoms(iatom2);
                if (atom2.element != constants::symbols::carbon) {continue;}

                double dist = atom1.distance(atom2);
                if (dist > min_dist) {continue;}

                min_dist = dist;
                min_atom1 = iatom1;
                min_atom2 = iatom2;
            }
        }

        constraints.emplace_back(std::make_shared<DistanceConstraint>(manager->protein, ibody1, ibody2, min_atom1, min_atom2));
        if (setting::general::verbose) {
            std::cout << "\tConstraint created between bodies " << ibody1 << " and " << ibody2 << " on atoms " << body1.atoms(min_atom1).name << " and " << body2.atoms(min_atom2).name << std::endl;
        }
    }

    if (constraints.empty()) {
        throw except::unexpected("rigidbody::constraints::generation::LinearConstraints::generate: No constraints were generated. This is probably a bug.");
    }

    return constraints;
}