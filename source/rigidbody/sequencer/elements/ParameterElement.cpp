// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/ParameterElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/parameters/decay/DecayFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/Rigidbody.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

ParameterElement::ParameterElement(observer_ptr<LoopElement> owner, std::unique_ptr<parameter::ParameterGenerationStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}

ParameterElement::~ParameterElement() = default;

void ParameterElement::run() {
    owner->_get_rigidbody()->parameter_generator = strategy;
}

ParameterElement& ParameterElement::decay_strategy(std::unique_ptr<parameter::decay::DecayStrategy> strategy) {
    this->strategy->set_decay_strategy(std::move(strategy));
    return *this;
}

ParameterElement& ParameterElement::max_rotation_angle(double radians) {
    strategy->set_max_rotation_angle(radians);
    return *this;
}

ParameterElement& ParameterElement::max_translation_distance(double distance) {
    strategy->set_max_translation_distance(distance);
    return *this;
}

observer_ptr<parameter::ParameterGenerationStrategy> ParameterElement::get_parameter_strategy() const {
    return strategy.get();
}

namespace {
    struct ParameterStrategyDefs {
        static inline std::string ROTATE_ONLY = "rotate_only";
        static inline std::string TRANSLATE_ONLY = "translate_only";
        static inline std::string BOTH = "both";
        static inline std::string SYMMETRY_ONLY = "symmetry_only";
    };
}

std::unique_ptr<GenericElement> ParameterElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    enum class Args {iterations, translate, rotate, strategy, decay_strategy};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::iterations,      {"iterations"}},
        {Args::translate,       {"translate"}},
        {Args::rotate,          {"rotate"}},
        {Args::strategy,        {"mode", "strategy"}},
        {Args::decay_strategy,  {"decay"}}
    };

    static auto get_parameter_strategy = [] (std::string_view line) {
        if (line == ParameterStrategyDefs::ROTATE_ONLY) {return settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly;}
        if (line == ParameterStrategyDefs::TRANSLATE_ONLY) {return settings::rigidbody::ParameterGenerationStrategyChoice::TranslationsOnly;}
        if (line == ParameterStrategyDefs::SYMMETRY_ONLY) {return settings::rigidbody::ParameterGenerationStrategyChoice::SymmetryOnly;}
        if (line == ParameterStrategyDefs::BOTH || line == "rotate_and_translate") {return settings::rigidbody::ParameterGenerationStrategyChoice::Simple;}
        throw except::parse_error("parameter", "Unknown choice \"" + std::string(line) + "\"");
    };

    static auto get_decay_strategy = [] (std::string_view line) {
        if (line == "linear") {return settings::rigidbody::DecayStrategyChoice::Linear;}
        if (line == "exponential") {return settings::rigidbody::DecayStrategyChoice::Exponential;}
        if (line == "none") {return settings::rigidbody::DecayStrategyChoice::None;}
        throw except::parse_error("parameter", "Unknown choice \"" + std::string(line) + "\"");
    };

    auto iterations = args.get<int>(valid_args[Args::iterations]);
    auto translate = args.get<double>(valid_args[Args::translate], 0);
    auto rotate = args.get<double>(valid_args[Args::rotate], 0);
    auto strategy = args.get<std::string>(valid_args[Args::strategy], "both");
    auto decay_strategy = args.get<std::string>(valid_args[Args::decay_strategy], "linear");

    bool has_translate = translate.found && translate.value != 0;
    bool has_rotate = rotate.found && rotate.value != 0;

    // validate arguments
    if (!args.inlined.empty()) {throw except::parse_error("parameter", "Unexpected inline arguments.");}
    if (!iterations.found) {throw except::parse_error("parameter", "Missing required argument \"iterations\".");}
    if (!strategy.found) {
        // deduce strategy from presence/absence of args

        // translate != 0; rotate = 0  --> translate_only
        if (has_translate && !has_rotate) {strategy.value = ParameterStrategyDefs::TRANSLATE_ONLY;}

        // translate = 0; rotate != 0  --> rotate_only
        else if (!has_translate && has_rotate) {strategy.value = ParameterStrategyDefs::ROTATE_ONLY;}

        // translate == 0; rotate = 0 --> symmetry_only (if symmetries are present)
        else if (!has_translate && !has_rotate) {
            if (owner->_get_rigidbody()->molecule.symmetry().has_symmetries()) {
                strategy.value = ParameterStrategyDefs::SYMMETRY_ONLY;
            } else {
                throw except::parse_error("parameter", "Missing one of \"mode\", \"translate\", or \"rotate\".");
            }
        }
    } else { // explicit strategy provided, so validate that required args are present and that no extraneous args are provided
        if (strategy.value == ParameterStrategyDefs::TRANSLATE_ONLY && !has_translate) {
            throw except::parse_error("parameter", "Missing required argument \"translate\" for strategy \"translate_only\".");
        } else if (strategy.value == ParameterStrategyDefs::ROTATE_ONLY && !has_rotate) {
            throw except::parse_error("parameter", "Missing required argument \"rotate\" for strategy \"rotate_only\".");
        } else if (strategy.value == ParameterStrategyDefs::BOTH && !(has_translate && has_rotate)) {
            throw except::parse_error("parameter", "Missing required arguments \"translate\" and \"rotate\" for strategy \"both\".");
        } else if (strategy.value == ParameterStrategyDefs::SYMMETRY_ONLY && (has_translate || has_rotate)) {
            throw except::parse_error("parameter", "Unexpected arguments \"translate\" and/or \"rotate\" for strategy \"symmetry_only\".");
        }
    }

    return std::make_unique<ParameterElement>(
        owner,
        rigidbody::factory::create_parameter_strategy(
            owner->_get_rigidbody(),
            rigidbody::factory::create_decay_strategy(iterations.value, get_decay_strategy(decay_strategy.value)),
            translate.value,
            rotate.value,
            get_parameter_strategy(strategy.value)
        )
    );
}