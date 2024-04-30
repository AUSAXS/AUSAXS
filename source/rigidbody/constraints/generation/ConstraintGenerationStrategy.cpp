/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

using namespace rigidbody::constraints;

ConstraintGenerationStrategy::ConstraintGenerationStrategy(const ConstraintManager* manager) : manager(manager) {}
ConstraintGenerationStrategy::~ConstraintGenerationStrategy() = default;