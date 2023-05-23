#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/RigidBody.h>
#include <utility/Exceptions.h>

#include <utility>

using namespace rigidbody::selection;

RandomConstraintSelect::RandomConstraintSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {
    unsigned int M = rigidbody->constraints->distance_constraints.size();
    std::random_device random;
    generator = std::mt19937(random());
    distribution = std::uniform_int_distribution<int>(0, M-1);
}

RandomConstraintSelect::~RandomConstraintSelect() = default;

std::pair<unsigned int, unsigned int> RandomConstraintSelect::next() {
    unsigned int iconstraint = distribution(generator);
    auto constraint = rigidbody->constraints->distance_constraints[iconstraint];
    unsigned int ibody = constraint.ibody1;
    for (unsigned int i = 0; i < rigidbody->constraints->distance_constraints_map.at(ibody).size(); i++) {
        if (rigidbody->constraints->distance_constraints_map.at(ibody).at(i) == &constraint) {
            return std::make_pair(ibody, i);
        }
    }
    throw except::invalid_argument("RandomConstraintSelect::next: Constraint " + std::to_string(iconstraint) + " not found");
}