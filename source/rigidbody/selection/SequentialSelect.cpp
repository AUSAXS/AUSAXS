#include <rigidbody/selection/SequentialSelect.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody;

SequentialSelect::SequentialSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {
    std::random_device random;
    generator = std::mt19937(random());
    distribution = std::uniform_int_distribution<int>(0, N-1);
}

SequentialSelect::~SequentialSelect() = default;

unsigned int SequentialSelect::next() {
    return distribution(generator);
}
