// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/AdditionalElements.h>
#include <rigidbody/sequencer/detail/ArgumentHelper.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/MessageElement.h>
#include <utility/StringUtils.h>
#include <utility/Random.h>

using namespace ausaxs::rigidbody::sequencer;

namespace {
    enum class OverlapArgs {scaling, distance};
    std::unordered_map<OverlapArgs, std::vector<std::string>> overlap_args_map = {
        {OverlapArgs::scaling, {"scaling", "factor"}},
        {OverlapArgs::distance, {"max", "max_distance", "distance"}}
    };
}

// seed [seed] - integer seed for the global random number generator
void detail::SeedElement::_parse(observer_ptr<LoopElement>, ParsedArgs&& args) {
    if (!utility::isinteger(args.inlined[0])) {throw except::parse_error("seed", "Expected an integer seed value, but got \"" + args.inlined[0] + "\".");}

    int seed = std::stoi(args.inlined[0]);
    random::set_seed(seed);
}

// log [message] - like print, but written to the log instead of the console
std::unique_ptr<GenericElement> detail::LogElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    auto message = args.inlined[0];
    return std::make_unique<MessageElement>(owner->_get_sequencer(), message, true);
}

// end - closes the innermost open block
void detail::LoopEndElement::_parse(observer_ptr<LoopElement>, ParsedArgs&&) {
}

// overlap_strength { scaling [factor], max [distance] }
void detail::OverlapStrengthElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    auto scaling = args.get<double>(overlap_args_map[OverlapArgs::scaling]);
    auto distance = args.get<double>(overlap_args_map[OverlapArgs::distance]);

    // inlined variant: overlap_strength strength distance
    if (!args.inlined.empty()) {
        if (args.inlined.size() == 1) {
            scaling = {std::stod(args.inlined[0]), true};
        } if (args.inlined.size() == 2) {
            distance = {std::stod(args.inlined[1]), true};
        }
    } else {
        if (!scaling.found) {throw except::parse_error("overlap_strength", "Missing required argument \"scaling\".");}
        if (!distance.found) {throw except::parse_error("overlap_strength", "Missing required argument \"distance\".");}
    }

    owner->_get_sequencer()->setup().set_overlap_function([a=scaling.value, d=distance.value] (double x) {return x < d ? a*std::pow((d-x)/d, 2) : 0;});
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
InlineSignature detail::OverlapStrengthElement::_valid_inline_arguments() { return {.names = {"strength", "distance"}, .min = 1, .max = 2}; }