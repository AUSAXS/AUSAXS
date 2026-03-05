// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/BodySelectElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::sequencer;

BodySelectElement::BodySelectElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}
BodySelectElement::~BodySelectElement() = default;

void BodySelectElement::run() {
    owner->_get_rigidbody()->body_selector = strategy;
}

std::unique_ptr<GenericElement> BodySelectElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    enum class Args {strategy, parameters};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::strategy,   {"point", "body"}},
        {Args::parameters, {"parameters", "parameter_mask", "mask"}},
    };
    auto strategy = args.get<std::string>(valid_args[Args::strategy]);
    auto mask_arg = args.get<std::string>(valid_args[Args::parameters]);

    if (!args.inlined.empty()) {throw except::parse_error("select", "Unexpected inline argument.");}
    if (args.named.size() < 1 || 2 < args.named.size()) {
        throw except::parse_error("select", "Invalid number of arguments. Expected 1 or 2, but got " + std::to_string(args.named.size()) + ".");
    }

    static auto get_body_select_strategy = [] (std::string_view line) {
        if (line == "random_body") {return settings::rigidbody::BodySelectStrategyChoice::RandomBodySelect;}
        if (line == "random_constraint") {return settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect;}
        if (line == "sequential_body") {return settings::rigidbody::BodySelectStrategyChoice::SequentialBodySelect;}
        if (line == "sequential_constraint") {return settings::rigidbody::BodySelectStrategyChoice::SequentialConstraintSelect;}
        throw except::parse_error("select", "Unknown choice \"" + std::string(line) + "\"");
    };

    static auto get_parameter_mask_strategy = [] (std::string_view line) {
        if (line == "all")                  {return settings::rigidbody::ParameterMaskStrategyChoice::All;}
        if (line == "real")                 {return settings::rigidbody::ParameterMaskStrategyChoice::Real;}
        if (line == "symmetry")             {return settings::rigidbody::ParameterMaskStrategyChoice::Symmetry;}
        if (line == "sequential")           {return settings::rigidbody::ParameterMaskStrategyChoice::Sequential;}
        if (line == "sequential_symmetry")  {return settings::rigidbody::ParameterMaskStrategyChoice::SequentialSymmetry;}
        if (line == "sequential_real")      {return settings::rigidbody::ParameterMaskStrategyChoice::SequentialReal;}
        if (line == "random")               {return settings::rigidbody::ParameterMaskStrategyChoice::Random;}
        throw except::parse_error("select", "Unknown mask strategy \"" + std::string(line) + "\"");
    };

    settings::rigidbody::BodySelectStrategyChoice body_strategy = strategy.found ? get_body_select_strategy(strategy.value) : settings::rigidbody::body_select_strategy;
    settings::rigidbody::ParameterMaskStrategyChoice mask_strategy = mask_arg.found ? get_parameter_mask_strategy(mask_arg.value) : settings::rigidbody::parameter_mask_strategy;
    return std::make_unique<BodySelectElement>(
        owner,
        rigidbody::factory::create_selection_strategy(
            owner->_get_rigidbody(),
            body_strategy,
            mask_strategy
        )
    );
}