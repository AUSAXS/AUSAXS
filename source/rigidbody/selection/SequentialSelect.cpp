#include <rigidbody/selection/SequentialSelect.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody;

SequentialSelect::SequentialSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialSelect::~SequentialSelect() = default;

std::pair<unsigned int, unsigned int> SequentialSelect::next() {
    unsigned int N = rigidbody->constraints->distance_constraints_map.size();
    unsigned int M = rigidbody->constraints->distance_constraints_map.at(ibody).size();

    if (iconstraint == M) {
        ibody = (ibody + 1) % N;
        iconstraint = 0;
    }

    return std::make_pair(ibody, iconstraint++);
}
