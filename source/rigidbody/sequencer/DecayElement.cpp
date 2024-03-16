/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/DecayElement.h>
#include <rigidbody/sequencer/ParameterElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

DecayElement::DecayElement(ParameterElement* owner) : ParameterElementCallback(owner), strategy(settings::rigidbody::decay_strategy) {}
DecayElement::DecayElement(ParameterElement* owner, const settings::rigidbody::DecayStrategyChoice& strategy) : ParameterElementCallback(owner), strategy(strategy) {}
DecayElement::~DecayElement() = default;

void DecayElement::apply() {
    std::cout << "DecayElement::apply()" << std::endl;
}