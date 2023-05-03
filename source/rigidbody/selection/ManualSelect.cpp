#include <rigidbody/selection/ManualSelect.h>

using namespace rigidbody::selection;

ManualSelect::ManualSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {}

ManualSelect::~ManualSelect() = default;

std::pair<unsigned int, unsigned int> ManualSelect::next() {}
