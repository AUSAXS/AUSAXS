#include <rigidbody/selection/ManualSelect.h>

using namespace rigidbody::selection;

ManualSelect::ManualSelect(const RigidBody* rigidbody);

ManualSelect::~ManualSelect() override;

std::pair<unsigned int, unsigned int> ManualSelect::next() override;
