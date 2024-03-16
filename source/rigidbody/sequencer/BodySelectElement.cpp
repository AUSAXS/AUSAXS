/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/BodySelectElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

BodySelectElement::BodySelectElement(LoopElement* owner) : LoopElementCallback(owner), strategy(settings::rigidbody::body_select_strategy) {}
BodySelectElement::BodySelectElement(LoopElement* owner, settings::rigidbody::BodySelectStrategyChoice strategy) : LoopElementCallback(owner), strategy(strategy) {}
BodySelectElement::~BodySelectElement() = default;

void BodySelectElement::apply() {
    std::cout << "BodySelectElement::apply()" << std::endl;
}