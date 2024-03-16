/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/selection/SequentialSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody::selection;

SequentialSelect::SequentialSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialSelect::~SequentialSelect() = default;

std::pair<unsigned int, unsigned int> SequentialSelect::next() {
    unsigned int N = rigidbody->get_constraint_manager()->distance_constraints_map.size();
    unsigned int M = rigidbody->get_constraint_manager()->distance_constraints_map.at(ibody).size();

    if (iconstraint == M) {
        ibody = (ibody + 1) % N;
        iconstraint = 0;
    }

    return std::make_pair(ibody, iconstraint++);
}
