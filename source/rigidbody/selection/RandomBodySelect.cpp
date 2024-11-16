/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/RigidBody.h>
#include <utility/Exceptions.h>

using namespace ausaxs::rigidbody::selection;

RandomBodySelect::RandomBodySelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {
    std::random_device random;
    generator = std::mt19937(random());
    distribution = std::uniform_int_distribution<int>(0, N-1);
}

RandomBodySelect::~RandomBodySelect() = default;

std::pair<unsigned int, int> RandomBodySelect::next() {
    unsigned int ibody = distribution(generator);

    unsigned int N = rigidbody->get_constraint_manager()->distance_constraints_map.at(ibody).size();
    switch (N) {
        case 0: {
            return std::make_pair(ibody, -1);
        }
        case 1: {
            return std::make_pair(ibody, 0);
        }
        default: {
            std::random_device random;
            std::mt19937 generator2(random());
            std::uniform_int_distribution<int> distribution2(0, N-1);
            unsigned int iconstraint = distribution2(generator2);

            return std::make_pair(ibody, iconstraint);
        }
    }
}