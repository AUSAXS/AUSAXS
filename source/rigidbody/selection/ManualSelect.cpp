#include <rigidbody/selection/ManualSelect.h>
#include <utility/Exceptions.h>

using namespace rigidbody::selection;

ManualSelect::ManualSelect(const RigidBody* rigidbody) : BodySelectStrategy(rigidbody) {}

ManualSelect::~ManualSelect() = default;

std::pair<unsigned int, unsigned int> ManualSelect::next() {
    throw except::not_implemented("ManualSelect::next: Not implemented.");
}
