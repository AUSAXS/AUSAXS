/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/selection/BodySelectStrategy.h>

#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::selection;

BodySelectStrategy::BodySelectStrategy(observer_ptr<const Rigidbody> rigidbody) : rigidbody(rigidbody), N(rigidbody->molecule.size_body()) {}