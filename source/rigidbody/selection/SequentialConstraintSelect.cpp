/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/RigidBody.h>

using namespace ausaxs::rigidbody::selection;

SequentialConstraintSelect::SequentialConstraintSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialConstraintSelect::~SequentialConstraintSelect() = default;

std::pair<unsigned int, int> SequentialConstraintSelect::next() {
    unsigned int N = rigidbody->get_constraint_manager()->distance_constraints_map.size();
    unsigned int M = rigidbody->get_constraint_manager()->distance_constraints_map.at(ibody).size();

    if (iconstraint == M) {
        ibody = (ibody + 1) % N;
        iconstraint = 0;
    }

    return std::make_pair(ibody, iconstraint++);
}
