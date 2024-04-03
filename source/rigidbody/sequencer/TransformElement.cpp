/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/TransformElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

TransformElement::TransformElement(LoopElement* owner) : LoopElementCallback(owner), strategy(settings::rigidbody::transform_strategy) {}
TransformElement::TransformElement(LoopElement* owner, settings::rigidbody::TransformationStrategyChoice strategy) : LoopElementCallback(owner), strategy(strategy) {}
TransformElement::~TransformElement() = default;

void TransformElement::run() {
    std::cout << "TransformElement::apply()" << std::endl;
}