/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/TransformElement.h>
#include <rigidbody/sequencer/RigidBodyManager.h>

#include <iostream>

using namespace rigidbody::sequencer;

TransformElement::TransformElement(LoopElement* owner, std::unique_ptr<rigidbody::transform::TransformStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}
TransformElement::~TransformElement() = default;

void TransformElement::run() {
    std::cout << "TransformElement::apply()" << std::endl;
    rigidbody->set_transform_manager(strategy);
}