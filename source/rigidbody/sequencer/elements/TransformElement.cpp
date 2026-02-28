// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/TransformElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/Rigidbody.h>
#include <settings/RigidBodySettings.h>

using namespace ausaxs::rigidbody::sequencer;

TransformElement::TransformElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::transform::TransformStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}
TransformElement::~TransformElement() = default;

void TransformElement::run() {
    owner->_get_rigidbody()->transformer = strategy;
}

std::unique_ptr<GenericElement> TransformElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("transform", "Unexpected named arguments.");}
    if (args.inlined.size() != 1) {throw except::parse_error("transform", "Expected only a single inline argument.");}

    static auto get_transform_strategy = [] (std::string_view line) {
        if (line == "rigid_transform" || line == "rigid") {return settings::rigidbody::TransformationStrategyChoice::RigidTransform;}
        if (line == "single_transform" || line == "single") {return settings::rigidbody::TransformationStrategyChoice::SingleTransform;}
        throw except::parse_error("transform", "Unknown choice \"" + std::string(line) + "\"");
    };

    settings::rigidbody::TransformationStrategyChoice strategy = get_transform_strategy(args.inlined[0]);
    return std::make_unique<TransformElement>(
        owner,
        rigidbody::factory::create_transform_strategy(
            owner->_get_rigidbody(),
            strategy
        )
    );
}