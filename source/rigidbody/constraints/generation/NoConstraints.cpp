/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/constraints/generation/NoConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>

using namespace rigidbody::constraints;

std::vector<DistanceConstraint> NoConstraints::generate() const {return {};}