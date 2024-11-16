/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/setup/RelativeHydrationElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/RigidBody.h>
#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/culling/CullingFactory.h>
#include <hydrate/culling/BodyCounterCulling.h>
#include <hydrate/generation/GridBasedHydration.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

RelativeHydrationElement::RelativeHydrationElement(observer_ptr<Sequencer> owner, const std::vector<double>& ratios) : owner(owner), ratios(ratios) {}
RelativeHydrationElement::RelativeHydrationElement(observer_ptr<Sequencer> owner, const std::vector<double>& ratios, const std::vector<std::string>& names) : owner(owner), ratios(ratios) {
    if (names.size() != ratios.size()) {
        throw std::runtime_error("RelativeHydrationElement::RelativeHydrationElement: The number of names and ratios must be equal.");
    }

    for (unsigned int i = 0; i < names.size(); ++i) {
        if (!owner->_get_body_names().contains(names[i])) {
            throw std::runtime_error("RelativeHydrationElement::RelativeHydrationElement: The body name \"" + names[i] + "\" is not known.");
        }
        this->ratios.push_back(owner->_get_body_names().at(names[i]));
    }
}

RelativeHydrationElement::~RelativeHydrationElement() = default;

void RelativeHydrationElement::run() {
    auto culling_strategy = hydrate::factory::construct_culling_strategy(owner->_get_rigidbody(), settings::hydrate::CullingStrategy::RandomBodyCounterStrategy);
    static_cast<hydrate::BodyCounterCulling*>(culling_strategy.get())->set_body_ratios(ratios);

    assert(
        dynamic_cast<hydrate::GridBasedHydration*>(owner->_get_rigidbody()->get_hydration_generator()) != nullptr && 
        "RelativeHydrationElement::run: owner->_get_rigidbody()->get_hydration_generator() is not a GridBasedHydration"
    );

    static_cast<hydrate::GridBasedHydration*>(owner->_get_rigidbody()->get_hydration_generator())->set_culling_strategy(std::move(culling_strategy));
    owner->_get_rigidbody()->generate_new_hydration();
}