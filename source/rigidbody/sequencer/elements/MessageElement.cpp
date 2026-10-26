// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/MessageElement.h>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/ArgumentHelper.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <utility/Console.h>
#include <utility/Logging.h>
#include <utility/StringUtils.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

std::function<std::string()> MessageElement::parse_user_msg(std::string_view msg) const {
    static auto iteration = [] () { return std::to_string(LoopElement::_get_current_iteration()); };
    static auto iterations_total  = [] () { return std::to_string(LoopElement::_get_total_iterations()); };
    auto dof = [this] () { return static_cast<double>(owner->_get_rigidbody()->controller->get_fitter()->dof()); };
    auto chi2 = [this, dof] () { return utility::round_double(owner->_get_current_conf()->chi2/dof(), 3); };
    auto chi2_best = [this, dof] () { return utility::round_double(owner->_get_best_conf()->chi2/dof(), 3); };
    auto chi2_penalty = [this, dof] () { return utility::round_double(owner->_get_rigidbody()->controller->get_fitter()->constraints_chi2()/dof(), 3); };
    auto chi2_no_penalty = [this, dof] () { return utility::round_double((owner->_get_current_conf()->chi2 - owner->_get_rigidbody()->controller->get_fitter()->constraints_chi2())/dof(), 3);};
    auto body_names = [this] () {
        const auto& entries = owner->_get_sequencer()->_body_name_registry().all();

        std::vector<std::string> parts;
        parts.reserve(entries.size());
        for (const auto& [index, entry] : entries) { // in index order, so each body is followed by its own symmetry replicas
            parts.push_back(entry.alias.empty() ? entry.default_name : entry.alias + " (" + entry.default_name + ")");
        }
        return utility::join(parts, ", ");
    };

    auto known_placeholder = [&] (std::string_view placeholder) -> std::function<std::string()> {
        if (placeholder == "{iteration}") {return iteration;}
        if (placeholder == "{iterations_total}") {return iterations_total;}
        if (placeholder == "{chi2}") {return chi2;}
        if (placeholder == "{chi2_best}") {return chi2_best;}
        if (placeholder == "{chi2_penalty}") {return chi2_penalty;}
        if (placeholder == "{chi2_no_penalty}") {return chi2_no_penalty;}
        if (placeholder == "{body_names}") {return body_names;}
        return nullptr;
    };

    // parts[i] is the text before funcs[i], and parts always has one more entry than funcs. A new part is started only
    // at a known placeholder: an unknown one, or a lone '{', is ordinary text and stays in the current part.
    std::vector<std::string> parts(1);
    std::vector<std::function<std::string()>> funcs;
    for (std::size_t i = 0; i < msg.size(); ++i) {
        if (msg[i] == '{') {
            if (std::size_t end = msg.find('}', i); end != std::string_view::npos) {
                if (auto func = known_placeholder(msg.substr(i, end - i + 1))) {
                    funcs.emplace_back(std::move(func));
                    parts.emplace_back();
                    i = end;
                    continue;
                }
            }
        }
        parts.back() += msg[i];
    }
    assert(parts.size() == funcs.size() + 1);
    return [parts, funcs] () {
        std::string result = parts[0];
        for (int i = 0; i < static_cast<int>(funcs.size()); ++i) {
            result += funcs[i]();
            result += parts[i+1];
        }
        return result;
    };
}

MessageElement::MessageElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, std::string_view message, bool log) 
    : MessageElement(owner, message, "white", log) 
{}

MessageElement::MessageElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, std::string_view message, std::string_view colour, bool log) 
    : owner(owner)
{
    message_func = [log, colour=console::color::parse(colour), builder=parse_user_msg(message)] () -> void {
        if (log) {
            logging::log(builder());
        } else {
            console::print_text(builder(), colour);
        }
    };
}

MessageElement::~MessageElement() = default;

void MessageElement::run() {
    assert(message_func && "MessageElement::run: message_func is not set.");
    message_func();
}

namespace {
    enum class Args {message, color};
    std::unordered_map<Args, std::vector<std::string>> args_map = {
        {Args::message, {"message", "msg"}},
        {Args::color, {"color", "colour"}}
    };
}

std::vector<std::string> MessageElement::_valid_arguments() {
    static auto map = detail::get_arg_names<Args>(args_map);
    return map;
}

InlineSignature MessageElement::_valid_inline_arguments() {
    return {.names = {"color", "message"}, .min = 0, .max = 2};
}

std::unique_ptr<GenericElement> MessageElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) { // NOLINT
    auto message = args.get<std::string>(args_map[Args::message]);
    auto color = args.get<std::string>(args_map[Args::color], "white");

    // inlined usage patterns
    if (!args.inlined.empty()) {
        // pattern 1: [message]
        if (args.inlined.size() == 1) {
            message.value = args.inlined[0];
            message.found = true;
        } else { // pattern 2: [color] [message]
            color.value = args.inlined[0];
            color.found = true;
            message.value = args.inlined[1];
            message.found = true;
        }
    }
    if (!message.found) {throw except::parse_error("print", "Missing required argument \"message\".");}
    try {console::color::parse(color.value);}
    catch (const std::exception&) {throw except::parse_error("print", "Unknown colour \"" + color.value + "\".");}
    return std::make_unique<MessageElement>(owner->_get_sequencer(), message.value, color.value, false);
}