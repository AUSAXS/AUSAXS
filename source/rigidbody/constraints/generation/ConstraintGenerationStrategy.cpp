#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

rigidbody::ConstraintGenerationStrategy::ConstraintGenerationStrategy(const ConstraintManager* manager) : manager(manager) {}
rigidbody::ConstraintGenerationStrategy::~ConstraintGenerationStrategy() = default;