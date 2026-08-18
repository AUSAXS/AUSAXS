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
    struct ParameterStrategyDefs {
        static inline std::string ROTATE_ONLY = "rotate_only";
        static inline std::string TRANSLATE_ONLY = "translate_only";
        static inline std::string BOTH = "both";
        static inline std::string SYMMETRY = "symmetry";
    };
}

namespace {
    enum class Args {iterations, translate, rotate, sym_translate, sym_rotate, strategy, decay_strategy};
    std::unordered_map<Args, std::vector<std::string>> args_map = {
        {Args::iterations,      {"iterations"}},
        {Args::translate,       {"translate"}},
        {Args::rotate,          {"rotate"}},
        {Args::sym_translate,   {"sym_translate"}},
        {Args::sym_rotate,      {"sym_rotate"}},
        {Args::strategy,        {"mode", "strategy"}},
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

// parameter { iterations [n], and any of: translate, rotate, sym_translate, sym_rotate, mode, decay }
std::unique_ptr<GenericElement> ParameterElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    static auto get_decay_strategy = [] (std::string_view line) {
        if (line == "linear") {return settings::rigidbody::DecayStrategyChoice::Linear;}
        if (line == "exponential") {return settings::rigidbody::DecayStrategyChoice::Exponential;}
        if (line == "none") {return settings::rigidbody::DecayStrategyChoice::None;}
        throw except::parse_error("parameter", "Unknown choice \"" + std::string(line) + "\"");
    };

    auto iterations = args.get<int>(args_map[Args::iterations]);
    auto translate = args.get<double>(args_map[Args::translate], 0);
    auto rotate = args.get<double>(args_map[Args::rotate], 0);
    auto sym_translate = args.get<double>(args_map[Args::sym_translate], 0);
    auto sym_rotate = args.get<double>(args_map[Args::sym_rotate], 0);
    auto mode = args.get<std::string>(args_map[Args::strategy]);
    auto decay_strategy = args.get<std::string>(args_map[Args::decay_strategy], "linear");

    if (!iterations.found) {throw except::parse_error("parameter", "Missing required argument \"iterations\".");}

    auto rigidbody = owner->_get_rigidbody();

    // a component is optimised exactly when its amplitude is non-zero, so an amplitude is all the generator needs
    parameter::ParameterAmplitudes amplitudes;

    // reject the amplitudes a mode has no use for rather than silently discarding them
    auto reject_unused = [] (const std::string& mode, std::initializer_list<std::pair<std::string, bool>> unused) {
        for (const auto& [name, found] : unused) {
            if (found) {
                throw except::parse_error("parameter", "Unexpected argument \"" + name + "\" for mode \"" + mode + "\".");
            }
        }
    };

    auto require_symmetries = [&rigidbody] () {
        if (!rigidbody->molecule.symmetry().has_symmetries()) {
            throw except::parse_error("parameter", "Symmetry optimisation was requested, but the molecule has no symmetries.");
        }
    };

    if (!mode.found) {
        // no mode given, so the amplitudes speak for themselves
        amplitudes = {
            .translation = translate.value,
            .rotation = rotate.value,
            .symmetry_translation = sym_translate.value,
            .symmetry_rotation = sym_rotate.value
        };

        // nothing at all was named: fall back to optimising the symmetries, which is the only thing that can be done without an amplitude
        if (amplitudes.translation == 0 && amplitudes.rotation == 0 && amplitudes.symmetry_translation == 0 && amplitudes.symmetry_rotation == 0) {
            if (!rigidbody->molecule.symmetry().has_symmetries()) {
                throw except::parse_error("parameter", "Missing one of \"mode\", \"translate\", \"rotate\", \"sym_translate\", or \"sym_rotate\".");
            }
            amplitudes.symmetry_translation = parameter::default_symmetry_translation(rigidbody);
            amplitudes.symmetry_rotation = parameter::default_symmetry_rotation();
        }
    }

    // "symmetry" retargets the two plain amplitudes onto the symmetry components, so that the common case of optimising only the symmetries
    // does not have to spell out the longer names. Naming them explicitly alongside it is therefore a contradiction.
    else if (mode.value == ParameterStrategyDefs::SYMMETRY || mode.value == "symmetry_only") {
        reject_unused(mode.value, {{"sym_translate", sym_translate.found}, {"sym_rotate", sym_rotate.found}});
        require_symmetries();
        if (!translate.found && !rotate.found) { // no amplitudes at all: optimise both components at their defaults
            amplitudes.symmetry_translation = parameter::default_symmetry_translation(rigidbody);
            amplitudes.symmetry_rotation = parameter::default_symmetry_rotation();
        } else {
            amplitudes.symmetry_translation = translate.value;
            amplitudes.symmetry_rotation = rotate.value;
        }
    }

    else if (mode.value == ParameterStrategyDefs::TRANSLATE_ONLY) {
        reject_unused(mode.value, {{"rotate", rotate.found}, {"sym_translate", sym_translate.found}, {"sym_rotate", sym_rotate.found}});
        if (translate.value == 0) {
            throw except::parse_error("parameter", "Missing required argument \"translate\" for mode \"translate_only\".");
        }
        amplitudes.translation = translate.value;
    }

    else if (mode.value == ParameterStrategyDefs::ROTATE_ONLY) {
        reject_unused(mode.value, {{"translate", translate.found}, {"sym_translate", sym_translate.found}, {"sym_rotate", sym_rotate.found}});
        if (rotate.value == 0) {
            throw except::parse_error("parameter", "Missing required argument \"rotate\" for mode \"rotate_only\".");
        }
        amplitudes.rotation = rotate.value;
    }

    else if (mode.value == ParameterStrategyDefs::BOTH || mode.value == "rotate_and_translate") {
        reject_unused(mode.value, {{"sym_translate", sym_translate.found}, {"sym_rotate", sym_rotate.found}});
        if (translate.value == 0 || rotate.value == 0) {
            throw except::parse_error("parameter", "Missing required arguments \"translate\" and \"rotate\" for mode \"both\".");
        }
        amplitudes.translation = translate.value;
        amplitudes.rotation = rotate.value;
    }

    else {throw except::parse_error("parameter", "Unknown choice \"" + mode.value + "\"");}

    if (amplitudes.symmetry_translation != 0 || amplitudes.symmetry_rotation != 0) {require_symmetries();}

    return std::make_unique<ParameterElement>(
        owner,
        rigidbody::factory::create_parameter_strategy(
            rigidbody,
            rigidbody::factory::create_decay_strategy(iterations.value, get_decay_strategy(decay_strategy.value)),
            amplitudes
        )
    );
}
