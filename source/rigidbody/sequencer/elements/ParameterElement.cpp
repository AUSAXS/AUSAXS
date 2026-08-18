// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/ParameterElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/detail/ArgumentHelper.h>
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

ParameterElement& ParameterElement::max_symmetry_rotation_angle(double radians) {
    strategy->set_max_symmetry_rotation_angle(radians);
    return *this;
}

ParameterElement& ParameterElement::max_symmetry_translation_distance(double distance) {
    strategy->set_max_symmetry_translation_distance(distance);
    return *this;
}

observer_ptr<parameter::ParameterGenerationStrategy> ParameterElement::get_parameter_strategy() const {
    return strategy.get();
}

namespace {
    enum class Args {iterations, translate, rotate, sym_translate, sym_rotate, decay_strategy};
    std::unordered_map<Args, std::vector<std::string>> args_map = {
        {Args::iterations,      {"iterations"}},
        {Args::translate,       {"translate"}},
        {Args::rotate,          {"rotate"}},
        {Args::sym_translate,   {"sym_translate"}},
        {Args::sym_rotate,      {"sym_rotate"}},
        {Args::decay_strategy,  {"decay"}}
    };
}

std::vector<std::string> ParameterElement::_valid_arguments() {
    static auto map = detail::get_arg_names<Args>(args_map);
    return map;
}

InlineSignature ParameterElement::_valid_inline_arguments() {
    return {.names = {}, .min = 0, .max = 0};
}

// parameter { iterations [n], and any of: translate, rotate, sym_translate, sym_rotate, decay }
std::unique_ptr<GenericElement> ParameterElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    static auto get_decay_strategy = [] (std::string_view line) {
        if (line == "linear") {return settings::rigidbody::DecayStrategyChoice::Linear;}
        if (line == "exponential") {return settings::rigidbody::DecayStrategyChoice::Exponential;}
        if (line == "none") {return settings::rigidbody::DecayStrategyChoice::None;}
        throw except::parse_error("parameter", "Unknown choice \"" + std::string(line) + "\"");
    };

    auto iterations = args.get<int>(args_map[Args::iterations]);
    auto decay_strategy = args.get<std::string>(args_map[Args::decay_strategy], "linear");

    if (!iterations.found) {throw except::parse_error("parameter", "Missing required argument \"iterations\".");}

    // a component is optimised exactly when its amplitude is non-zero, so the amplitudes are the whole configuration
    parameter::ParameterAmplitudes amplitudes = {
        .translation = args.get<double>(args_map[Args::translate], 0).value,
        .rotation = args.get<double>(args_map[Args::rotate], 0).value,
        .symmetry_translation = args.get<double>(args_map[Args::sym_translate], 0).value,
        .symmetry_rotation = args.get<double>(args_map[Args::sym_rotate], 0).value
    };

    auto rigidbody = owner->_get_rigidbody();
    bool has_symmetries = rigidbody->molecule.symmetry().has_symmetries();

    // nothing was named: fall back to optimising the symmetries, the only thing that can be done without an amplitude
    if (amplitudes.translation == 0 && amplitudes.rotation == 0 && amplitudes.symmetry_translation == 0 && amplitudes.symmetry_rotation == 0) {
        if (!has_symmetries) {
            throw except::parse_error("parameter", "Missing one of \"translate\", \"rotate\", \"sym_translate\", or \"sym_rotate\".");
        }
        amplitudes.symmetry_translation = parameter::default_symmetry_translation(rigidbody);
        amplitudes.symmetry_rotation = parameter::default_symmetry_rotation();
    }

    if ((amplitudes.symmetry_translation != 0 || amplitudes.symmetry_rotation != 0) && !has_symmetries) {
        throw except::parse_error("parameter", "Symmetry optimisation was requested, but the molecule has no symmetries.");
    }

    return std::make_unique<ParameterElement>(
        owner,
        rigidbody::factory::create_parameter_strategy(
            rigidbody,
            rigidbody::factory::create_decay_strategy(iterations.value, get_decay_strategy(decay_strategy.value)),
            amplitudes
        )
    );
}
