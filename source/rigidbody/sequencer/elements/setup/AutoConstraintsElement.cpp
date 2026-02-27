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

AutoConstraintsElement::AutoConstraintsElement(observer_ptr<Sequencer> owner, settings::rigidbody::ConstraintGenerationStrategyChoice strategy) : owner(owner), strategy(strategy) {}

void AutoConstraintsElement::run() {
    if (owner->_get_rigidbody() == nullptr) {throw std::runtime_error("AutoConstraintsElement::run: No body is currently loaded.");}
    owner->_get_rigidbody()->constraints->generate_constraints(rigidbody::factory::generate_constraints(owner->_get_rigidbody()->constraints.get(), strategy));
}

std::unique_ptr<GenericElement> AutoConstraintsElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    static auto get_constraint_strategy = [] (std::string_view line) {
        if (line == "none") {return settings::rigidbody::ConstraintGenerationStrategyChoice::None;}
        if (line == "linear") {return settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;}
        if (line == "volumetric") {return settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric;}
        throw except::parse_error("autoconstrain", "Unknown choice \"" + std::string(line) + "\"");
    };

    if (!args.named.empty()) {throw except::parse_error("autoconstrain", "Unexpected named argument \"" + args.named.begin()->first + "\".");}
    if (args.inlined.size() != 1) {throw except::parse_error("autoconstrain", "Invalid number of inline arguments. Expected 1, but got " + std::to_string(args.inlined.size()) + ".");}
    return std::make_unique<AutoConstraintsElement>(owner->_get_sequencer(), get_constraint_strategy(args.inlined[0]));
}