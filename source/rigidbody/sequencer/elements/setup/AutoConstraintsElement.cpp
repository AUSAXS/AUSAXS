// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/setup/AutoConstraintsElement.h>
#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <settings/RigidBodySettings.h>

using namespace ausaxs::rigidbody::sequencer;

AutoConstraintsElement::AutoConstraintsElement(observer_ptr<Sequencer> owner, settings::rigidbody::ConstraintGenerationStrategyChoice strategy) : owner(owner), strategy(strategy) {
    owner->_get_rigidbody()->constraints->generate_constraints(rigidbody::factory::generate_constraints(owner->_get_rigidbody()->constraints.get(), strategy));
}

void AutoConstraintsElement::run() {}

std::vector<std::string> AutoConstraintsElement::_valid_arguments() {
    return {};
}

InlineSignature AutoConstraintsElement::_valid_inline_arguments() {
    return {.names = {"strategy"}, .min = 1, .max = 1};
}

// autoconstrain [strategy] - one of: none, backbone
std::unique_ptr<GenericElement> AutoConstraintsElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    static auto get_constraint_strategy = [] (std::string_view line) {
        if (line == "none") {return settings::rigidbody::ConstraintGenerationStrategyChoice::None;}
        if (line == "backbone") {return settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone;}
        throw except::parse_error("autoconstrain", "Unknown choice \"" + std::string(line) + "\"");
    };

    return std::make_unique<AutoConstraintsElement>(owner->_get_sequencer(), get_constraint_strategy(args.inlined[0]));
}