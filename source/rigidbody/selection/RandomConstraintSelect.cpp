#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody;

RandomConstraintSelect::RandomConstraintSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {
    unsigned int M = rigidbody->get_constraints().size();
    std::random_device random;
    generator = std::mt19937(random());
    distribution = std::uniform_int_distribution<int>(0, M-1);
}

RandomConstraintSelect::~RandomConstraintSelect() = default;

std::pair<unsigned int, unsigned int> RandomConstraintSelect::next() {
    unsigned int iconstraint = distribution(generator);
    auto constraint = rigidbody->get_constraint(iconstraint);
    unsigned int ibody = constraint->ibody1;
    for (unsigned int i = 0; i < rigidbody->constraint_map.at(ibody).size(); i++) {
        if (rigidbody->constraint_map.at(ibody).at(i) == constraint) {
            return std::make_pair(ibody, i);
        }
    }
    throw except::invalid_argument("RandomConstraintSelect::next: Constraint " + std::to_string(iconstraint) + " not found");
}