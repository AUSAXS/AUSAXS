/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody::sequencer;

BodySelectElement::BodySelectElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}
BodySelectElement::~BodySelectElement() = default;

void BodySelectElement::run() {
    owner->_get_rigidbody()->set_body_select_manager(strategy);
}