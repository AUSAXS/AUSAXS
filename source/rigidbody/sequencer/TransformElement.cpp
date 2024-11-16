/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/TransformElement.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/RigidBody.h>

using namespace ausaxs::rigidbody::sequencer;

TransformElement::TransformElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::transform::TransformStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}
TransformElement::~TransformElement() = default;

void TransformElement::run() {
    owner->_get_rigidbody()->set_transform_manager(strategy);
}