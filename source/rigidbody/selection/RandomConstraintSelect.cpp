/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/RigidBody.h>
#include <utility/Exceptions.h>

#include <utility>

using namespace rigidbody::selection;

RandomConstraintSelect::RandomConstraintSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {
    unsigned int M = rigidbody->get_constraint_manager()->distance_constraints.size();
    std::random_device random;
    generator = std::mt19937(random());
    distribution = std::uniform_int_distribution<int>(0, M-1);
}

RandomConstraintSelect::~RandomConstraintSelect() = default;

std::pair<unsigned int, unsigned int> RandomConstraintSelect::next() {
    unsigned int iconstraint = distribution(generator);
    const auto& constraint = rigidbody->get_constraint_manager()->distance_constraints[iconstraint];
    unsigned int ibody = constraint.ibody1;

    // find the index of the constraint in the list of constraints for the body
    for (unsigned int i = 0; i < rigidbody->get_constraint_manager()->distance_constraints_map.at(ibody).size(); i++) {
        // address comparison since the DistanceConstraint comparison operator is a weak equality comparing only its contents
        if (&rigidbody->get_constraint_manager()->distance_constraints_map.at(ibody).at(i).get() == &constraint) {
            return std::make_pair(ibody, i);
        }
    }
    throw except::invalid_argument("RandomConstraintSelect::next: Constraint " + std::to_string(iconstraint) + " not found");
}