// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/AdditionalElements.h>

#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/sequencer/detail/ArgumentHelper.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/MessageElement.h>
#include <utility/Random.h>
#include <utility/StringUtils.h>

using namespace ausaxs::rigidbody::sequencer;

namespace {
    enum class OverlapArgs {scaling, distance};
    std::unordered_map<OverlapArgs, std::vector<std::string>> overlap_args_map = {
        {OverlapArgs::scaling, {"scaling", "factor"}},
        {OverlapArgs::distance, {"max", "max_distance", "distance"}}
    };
}

// seed [seed] - integer seed for the global random number generator
void detail::SeedElement::_parse(observer_ptr<LoopElement> /*unused*/, ParsedArgs&& args) { // NOLINT
    if (!utility::isinteger(args.inlined[0])) {throw except::parse_error("seed", "Expected an integer seed value, but got \"" + args.inlined[0] + "\".");}

    int seed = std::stoi(args.inlined[0]);
    random::set_seed(seed);
}

// log [message] - like print, but written to the log instead of the console
std::unique_ptr<GenericElement> detail::LogElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) { // NOLINT
    auto message = args.inlined[0];
    return std::make_unique<MessageElement>(owner->_get_sequencer(), message, true);
}

// end - closes the innermost open block
void detail::LoopEndElement::_parse(observer_ptr<LoopElement> /*unused*/, ParsedArgs&& /*unused*/) {}

// overlap_strength [scaling] [distance], or overlap_strength { scaling [factor], max [distance] }
//  with a distance, the weight is scaling*((distance - r)/distance)^2 below the distance and 0 beyond it
//  without one, the default weight exp(-5r) is multiplied by the scaling
void detail::OverlapStrengthElement::_parse(observer_ptr<LoopElement> /*unused*/, ParsedArgs&& args) { // NOLINT
    auto scaling = args.get<double>(overlap_args_map[OverlapArgs::scaling], 1);
    auto distance = args.get<double>(overlap_args_map[OverlapArgs::distance]);

    // inlined variant: overlap_strength scaling [distance]
    if (!args.inlined.empty()) {
        if (!utility::isnumeric(args.inlined[0])) {throw except::parse_error("overlap_strength", "Expected a number for the scaling, but got \"" + args.inlined[0] + "\".");}
        scaling = {.value=std::stod(args.inlined[0]), .found=true};
        if (args.inlined.size() == 2) {
            if (!utility::isnumeric(args.inlined[1])) {throw except::parse_error("overlap_strength", "Expected a number for the distance, but got \"" + args.inlined[1] + "\".");}
            distance = {.value=std::stod(args.inlined[1]), .found=true};
        }
    }
    if (!scaling.found && !distance.found) {throw except::parse_error("overlap_strength", "Expected a scaling, a distance, or both.");}

    if (!distance.found) {
        rigidbody::constraints::OverlapConstraint::set_overlap_function([a=scaling.value] (double x) {return a*std::exp(-5*x);});
        return;
    }
    rigidbody::constraints::OverlapConstraint::set_overlap_function([a=scaling.value, d=distance.value] (double x) {return x < d ? a*std::pow((d-x)/d, 2) : 0;});
}

std::vector<std::string> detail::SeedElement::_valid_arguments() { return {}; }
std::vector<std::string> detail::LoopEndElement::_valid_arguments() { return {}; }
std::vector<std::string> detail::LogElement::_valid_arguments() { return {}; }
std::vector<std::string> detail::OverlapStrengthElement::_valid_arguments() {
    static auto map = get_arg_names<OverlapArgs>(overlap_args_map);
    return map;
}

InlineSignature detail::SeedElement::_valid_inline_arguments() { return {.names = {"seed"}, .min = 1, .max = 1}; }
InlineSignature detail::LoopEndElement::_valid_inline_arguments() { return {}; }
InlineSignature detail::LogElement::_valid_inline_arguments() { return {.names = {"message"}, .min = 1, .max = 1}; }
InlineSignature detail::OverlapStrengthElement::_valid_inline_arguments() { return {.names = {"scaling", "distance"}, .min = 0, .max = 2}; }