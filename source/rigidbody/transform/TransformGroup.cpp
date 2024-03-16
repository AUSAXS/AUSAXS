/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/Body.h>

using namespace rigidbody::transform;

TransformGroup::TransformGroup(std::vector<data::Body*> bodies, std::vector<unsigned int> indices, const constraints::DistanceConstraint& target, Vector3<double> pivot) 
    : bodies(bodies), indices(indices), target(target), pivot(pivot) {}

TransformGroup::~TransformGroup() = default;