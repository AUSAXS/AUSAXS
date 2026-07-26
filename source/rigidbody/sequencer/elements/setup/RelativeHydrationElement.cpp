// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/RelativeHydrationElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/Rigidbody.h>
#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/culling/CullingFactory.h>
#include <hydrate/culling/BodyCounterCulling.h>
#include <hydrate/generation/GridBasedHydration.h>
#include <data/Molecule.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

namespace {
    enum class Options {Maximum, High, Normal, Low, Minimum};

    double to_value(Options opt) {
        switch (opt) {
            case Options::Maximum:  return 1.75;
            case Options::High:     return 1.5;
            case Options::Normal:   return 1.0;
            case Options::Low:      return 0.5;
            case Options::Minimum:  return 0.25;
        }
        return 0;
    }
}

RelativeHydrationElement::RelativeHydrationElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& names, const std::vector<double>& ratios) : owner(owner) {
    assert(names.size() == ratios.size() && "RelativeHydrationElement::RelativeHydrationElement: The number of names and ratios must be equal.");
    const auto& body_names = owner->setup()._body_name_registry();
    this->ratios.assign(owner->_get_molecule()->size_body(), to_value(Options::Normal));

    for (unsigned int i = 0; i < names.size(); ++i) {
        if (!body_names.contains(names[i])) {
            throw std::runtime_error("RelativeHydrationElement::RelativeHydrationElement: The body name \"" + names[i] + "\" is not known.");
        }
        int ibody = body_names.resolve_body(names[i]); // throws on a symmetry replica, which has no hydration of its own to scale
        assert(ibody < static_cast<int>(this->ratios.size()) && "RelativeHydrationElement::RelativeHydrationElement: body name registry is out of sync with the molecule.");
        this->ratios[static_cast<std::size_t>(ibody)] = ratios[i];
    }
}

RelativeHydrationElement::~RelativeHydrationElement() = default;

const std::vector<double>& RelativeHydrationElement::_get_ratios() const {
    return ratios;
}

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

std::vector<std::string> RelativeHydrationElement::_valid_arguments() {
    return {};
}

std::unique_ptr<GenericElement> RelativeHydrationElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    static const std::unordered_map<std::string, Options> options = {
        {"max",     Options::Maximum},
        {"maximum", Options::Maximum},
        {"high",    Options::High},
        {"normal",  Options::Normal},
        {"low",     Options::Low},
        {"minimum", Options::Minimum},
        {"min",     Options::Minimum}
    };

    // inline usage pattern: [body] [hydration level]
    const auto& body_names = owner->_get_sequencer()->setup()._body_name_registry();
    if (!args.inlined.empty()) {
        if (args.inlined.size() != 2) {throw except::parse_error("relative_hydration", "Only 2 inline arguments can be provided.");}
        if (!body_names.contains(args.inlined[0])) {throw except::parse_error("relative_hydration", "Unknown body name \"" + args.inlined[0] + "\".");}
        if (!options.contains(args.inlined[1])) {throw except::parse_error("relative_hydration", "Unknown hydration level \"" + args.inlined[1] + "\".");}
        return std::make_unique<RelativeHydrationElement>(
            owner->_get_sequencer(),
            std::vector<std::string>{args.inlined[0]},
            std::vector<double>{to_value(options.at(args.inlined[1]))}
        );
    }

    // named usage pattern: 
    // relative_hydration {
    //     [body] [hydration level]
    //     [body] [hydration level]
    //     ...
    // }
    std::vector<double> ratios;
    std::vector<std::string> names;
    for (const auto& [name, value] : args.named) {
        if (!body_names.contains(name)) {throw except::parse_error("relative_hydration", "Unknown body name \"" + name + "\".");}
        if (value.size() != 1) {throw except::parse_error("relative_hydration", "Expected only a single argument.");}
        if (!options.contains(value.args[0])) {throw except::parse_error("relative_hydration", "Unknown hydration level \"" + value.args[0].str + "\".");}
        names.emplace_back(name);
        ratios.push_back(to_value(options.at(value.args[0])));
    }
    return std::make_unique<RelativeHydrationElement>(owner->_get_sequencer(), names, ratios);
}