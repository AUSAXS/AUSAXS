// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/RelativeHydrationElement.h>
#include <rigidbody/Rigidbody.h>
#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/culling/CullingFactory.h>
#include <hydrate/culling/BodyCounterCulling.h>
#include <hydrate/generation/GridBasedHydration.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

RelativeHydrationElement::RelativeHydrationElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<double>& ratios) : owner(owner), ratios(ratios) {
    assert(names.size() == ratios.size() && "RelativeHydrationElement::RelativeHydrationElement: The number of names and ratios must be equal.");

    for (unsigned int i = 0; i < names.size(); ++i) {
        if (!owner->setup()._get_body_names().contains(names[i])) {
            throw std::runtime_error("RelativeHydrationElement::RelativeHydrationElement: The body name \"" + names[i] + "\" is not known.");
        }
        this->ratios.push_back(owner->setup()._get_body_names().at(names[i]));
    }
}

RelativeHydrationElement::~RelativeHydrationElement() = default;

void RelativeHydrationElement::run() {
    auto culling_strategy = hydrate::factory::construct_culling_strategy(owner->_get_molecule(), settings::hydrate::CullingStrategy::RandomCounterStrategy);
    static_cast<hydrate::BodyCounterCulling*>(culling_strategy.get())->set_body_ratios(ratios);

    assert(
        dynamic_cast<hydrate::GridBasedHydration*>(owner->_get_molecule()->get_hydration_generator()) != nullptr && 
        "RelativeHydrationElement::run: owner->_get_rigidbody()->get_hydration_generator() is not a GridBasedHydration"
    );

    static_cast<hydrate::GridBasedHydration*>(owner->_get_molecule()->get_hydration_generator())->set_culling_strategy(std::move(culling_strategy));
    owner->_get_molecule()->generate_new_hydration();
}