#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

using namespace rigidbody::constraints;

ConstraintGenerationStrategy::ConstraintGenerationStrategy(const ConstraintManager* manager) : manager(manager) {}
ConstraintGenerationStrategy::~ConstraintGenerationStrategy() = default;