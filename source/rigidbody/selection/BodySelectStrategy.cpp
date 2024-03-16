/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/selection/BodySelectStrategy.h>

#include <rigidbody/RigidBody.h>

using namespace rigidbody::selection;

BodySelectStrategy::BodySelectStrategy(const RigidBody* rigidbody) : rigidbody(rigidbody), N(rigidbody->body_size()) {}