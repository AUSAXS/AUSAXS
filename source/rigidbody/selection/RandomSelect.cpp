#include <rigidbody/selection/RandomSelect.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody;

RandomSelect::RandomSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {
    std::random_device random;
    generator = std::mt19937(random());
    distribution = std::uniform_int_distribution<int>(0, N-1);
}

RandomSelect::~RandomSelect() = default;

unsigned int RandomSelect::next() {
    return distribution(generator);
}

