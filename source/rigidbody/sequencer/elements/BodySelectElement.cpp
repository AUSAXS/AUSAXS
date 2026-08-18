// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/BodySelectElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/ArgumentHelper.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/Rigidbody.h>

#include <optional>

using namespace ausaxs::rigidbody::sequencer;

BodySelectElement::BodySelectElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}
BodySelectElement::~BodySelectElement() = default;

void BodySelectElement::run() {
    owner->_get_rigidbody()->body_selector = strategy;
}

namespace {
    enum class Args {strategy, parameters};
    std::unordered_map<Args, std::vector<std::string>> args_map = {
        {Args::strategy,   {"point", "body"}},
        {Args::parameters, {"parameters", "parameter_mask", "mask"}},
    };
}

std::vector<std::string> BodySelectElement::_valid_arguments() {
    static auto map = detail::get_arg_names<Args>(args_map);
    return map;
}

InlineSignature BodySelectElement::_valid_inline_arguments() {
    return {.names = {"strategy or body name"}, .min = 0, .max = 1};
}

namespace {
    std::optional<ausaxs::settings::rigidbody::BodySelectStrategyChoice> try_get_body_select_strategy(std::string_view line) {
        using Choice = ausaxs::settings::rigidbody::BodySelectStrategyChoice;
        if (line == "random_body") {return Choice::RandomBodySelect;}
        if (line == "random_constraint") {return Choice::RandomConstraintSelect;}
        if (line == "sequential_body") {return Choice::SequentialBodySelect;}
        if (line == "sequential_constraint") {return Choice::SequentialConstraintSelect;}
        return std::nullopt;
    }
}

std::unique_ptr<GenericElement> BodySelectElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    // inlined usage patterns
    if (!args.inlined.empty()) {
        const std::string& token = args.inlined[0];

        // pattern 1: [strategy]
        if (auto choice = try_get_body_select_strategy(token)) {
            return std::make_unique<BodySelectElement>(owner, rigidbody::factory::create_selection_strategy(owner->_get_rigidbody(), *choice));
        }

        const auto& body_names = owner->_get_sequencer()->setup()._body_name_registry();
        if (!body_names.contains(token)) {
            throw except::parse_error("select", "Unknown body select strategy, body name/alias, or symmetry tag \"" + token + "\".");
        }

        auto sel = body_names.resolve(token);

        // pattern 2: [body name or alias] - optimize its pose and all of its symmetries
        if (sel.symmetry == -1) {
            return std::make_unique<BodySelectElement>(owner, rigidbody::factory::create_manual_selection_strategy(owner->_get_rigidbody(), sel.body));
        }

        // pattern 3: [symmetry tag] (e.g. "b1s2") - optimize that one declared symmetry in isolation. The replica is deliberately
        // ignored, since all replicas of a symmetry are generated from its single shared parameter set - b1s2 and b1s2r3 target the
        // same thing. A tag naming a symmetry shared between several bodies likewise resolves to the body owning it, so any
        // participant's name reaches the same parameters.
        return std::make_unique<BodySelectElement>(
            owner,
            rigidbody::factory::create_manual_symmetry_selection_strategy(owner->_get_rigidbody(), sel.body, sel.symmetry)
        );
    }

    // named form: a select strategy and/or a parameter mask, each falling back to its global setting
    auto strategy = args.get<std::string>(args_map[Args::strategy]);
    auto mask_arg = args.get<std::string>(args_map[Args::parameters]);

    if (!strategy.found && !mask_arg.found) {
        throw except::parse_error("select", "Missing arguments. Expected a select strategy and/or a parameter mask.");
    }

    static auto get_body_select_strategy = [] (std::string_view line) {
        if (auto choice = try_get_body_select_strategy(line)) {return *choice;}
        throw except::parse_error("select", "Unknown choice \"" + std::string(line) + "\"");
    };

    static auto get_parameter_mask_strategy = [] (std::string_view line) {
        if (line == "all")                  {return settings::rigidbody::ParameterMaskStrategyChoice::All;}
        if (line == "real")                 {return settings::rigidbody::ParameterMaskStrategyChoice::Real;}
        if (line == "symmetry")             {return settings::rigidbody::ParameterMaskStrategyChoice::Symmetry;}
        if (line == "symmetry_translation") {return settings::rigidbody::ParameterMaskStrategyChoice::SymmetryTranslation;}
        if (line == "symmetry_axis")        {return settings::rigidbody::ParameterMaskStrategyChoice::SymmetryAxis;}
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