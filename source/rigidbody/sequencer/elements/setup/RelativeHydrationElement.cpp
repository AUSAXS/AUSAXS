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
#include <map>

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

    std::map<std::string, double> custom_levels;
    bool levels_installed = false;

    std::string permanent_name(const detail::BodyNameRegistry& registry, const std::string& name) {
        for (const auto& [index, entry] : registry.all()) {
            if (entry.default_name == name || entry.alias == name) {return entry.default_name;}
        }
        assert(false && "permanent_name: the registry holds no entry for a name it claims to know.");
        return name;
    }
}

RelativeHydrationElement::RelativeHydrationElement(observer_ptr<Sequencer> owner, const std::string& name, double ratio) : owner(owner) {
    const auto& body_names = owner->setup()._body_name_registry();
    if (!body_names.contains(name)) {
        throw std::runtime_error("RelativeHydrationElement::RelativeHydrationElement: The body name \"" + name + "\" is not known.");
    }
    static_cast<void>(body_names.resolve_body(name)); // throws on a symmetry replica, which has no hydration of its own to scale
    custom_levels[permanent_name(body_names, name)] = ratio;
}

RelativeHydrationElement::~RelativeHydrationElement() {
    custom_levels.clear();
    levels_installed = false;
}

std::vector<double> RelativeHydrationElement::_get_ratios() const {
    const auto& body_names = owner->setup()._body_name_registry();
    std::vector<double> ratios(owner->_get_molecule()->size_body(), to_value(Options::Normal));

    for (const auto& [name, ratio] : custom_levels) {
        if (!body_names.contains(name)) {
            throw std::runtime_error(
                "RelativeHydrationElement::_get_ratios: A relative hydration level was declared for body \"" + name + "\", but that body no longer exists."
            );
        }
        int ibody = body_names.resolve_body(name);
        assert(ibody < static_cast<int>(ratios.size()) && "RelativeHydrationElement::_get_ratios: body name registry is out of sync with the molecule.");
        ratios[static_cast<std::size_t>(ibody)] = ratio;
    }
    return ratios;
}

void RelativeHydrationElement::run() {
    if (levels_installed) {return;} // the first element to run installs every declaration; see custom_levels
    levels_installed = true;

    auto culling_strategy = hydrate::factory::construct_culling_strategy(owner->_get_molecule(), settings::hydrate::CullingStrategy::RandomCounterStrategy);
    static_cast<hydrate::BodyCounterCulling*>(culling_strategy.get())->set_body_ratios(_get_ratios());

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

    // usage pattern: [body] [hydration level]. To set a level on several bodies, repeat the element - the declarations
    // accumulate into one culling strategy, so nothing is lost and the hydration layer is still generated only once.
    if (args.inlined.size() != 2) {throw except::parse_error("relative_hydration", "Expected [body] [hydration level].");}

    const auto& body_names = owner->_get_sequencer()->setup()._body_name_registry();
    if (!body_names.contains(args.inlined[0])) {throw except::parse_error("relative_hydration", "Unknown body name \"" + args.inlined[0] + "\".");}
    if (!options.contains(args.inlined[1])) {throw except::parse_error("relative_hydration", "Unknown hydration level \"" + args.inlined[1] + "\".");}
    return std::make_unique<RelativeHydrationElement>(
        owner->_get_sequencer(),
        args.inlined[0],
        to_value(options.at(args.inlined[1]))
    );
}