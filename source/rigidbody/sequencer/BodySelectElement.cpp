/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/sequencer/RigidBodyManager.h>

#include <iostream>

using namespace rigidbody::sequencer;

BodySelectElement::BodySelectElement(LoopElement* owner, std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}
BodySelectElement::~BodySelectElement() = default;

void BodySelectElement::run() {
    std::cout << "BodySelectElement::apply()" << std::endl;
    rigidbody->set_body_select_manager(strategy);
}