// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/AdditionalElements.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/MessageElement.h>
#include <utility/StringUtils.h>
#include <utility/Random.h>

using namespace ausaxs::rigidbody::sequencer;

void detail::SeedElement::_parse(observer_ptr<LoopElement>, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("seed", "Unexpected named argument.");}
    if (args.inlined.size() != 1) {throw except::parse_error("seed", "Expected only a single inline argument.");}
    if (!utility::isinteger(args.inlined[0])) {throw except::parse_error("seed", "Expected an integer seed value, but got \"" + args.inlined[0] + "\".");}

    int seed = std::stoi(args.inlined[0]);
    random::set_seed(seed);
}

std::unique_ptr<GenericElement> detail::LogElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("log", "Unexpected named argument.");}
    if (args.inlined.size() != 1) {throw except::parse_error("log", "Expected only a single inline argument.");}
    auto message = args.inlined[0];
    return std::make_unique<MessageElement>(owner->_get_sequencer(), message, true);
}

void detail::LoopEndElement::_parse(observer_ptr<LoopElement>, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("end", "Unexpected named argument.");}
    if (!args.inlined.empty()) {throw except::parse_error("end", "Unexpected inline argument.");}
}

void detail::OverlapStrengthElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    enum class Args {scaling, distance};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::scaling, {"scaling", "factor"}},
        {Args::distance, {"max", "max_distance", "distance"}}
    };

    auto scaling = args.get<double>(valid_args[Args::scaling]);
    auto distance = args.get<double>(valid_args[Args::distance]);

    if (!scaling.found) {throw except::parse_error("overlap_strength", "Missing required argument \"scaling\".");}
    if (!distance.found) {throw except::parse_error("overlap_strength", "Missing required argument \"distance\".");}
    owner->_get_sequencer()->setup().set_overlap_function([a=scaling.value, d=distance.value] (double x) {return x < d ? a*std::pow((d-x)/d, 2) : 0;});
}