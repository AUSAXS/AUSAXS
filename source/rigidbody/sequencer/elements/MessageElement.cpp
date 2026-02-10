// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/MessageElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Console.h>
#include <utility/Logging.h>

using namespace ausaxs::rigidbody::sequencer;

std::function<std::string()> MessageElement::parse_user_msg(std::string_view msg) const {
    static auto iteration = [] () { return std::to_string(LoopElement::_get_current_iteration()); };
    static auto iterations_total  = [] () { return std::to_string(LoopElement::_get_total_iterations()); };
    auto chi2 = [this] () { return utility::round_double(owner->_get_current_conf()->chi2, 3); };
    auto chi2_best = [this] () { return utility::round_double(owner->_get_best_conf()->chi2, 3); };
    auto chi2_penalty = [this] () { return utility::round_double(owner->_get_rigidbody()->controller->get_fitter()->constraints_chi2(), 3); };
    auto chi2_no_penalty = [this] () { return utility::round_double(owner->_get_current_conf()->chi2 - owner->_get_rigidbody()->controller->get_fitter()->constraints_chi2(), 3);};

    std::vector<std::string> parts;
    std::vector<std::function<std::string()>> funcs;

    int last_pos = 0;
    for (size_t i = 0; i < msg.size(); ++i) {
        if (msg[i] == '{') {
            parts.emplace_back(msg.substr(last_pos, i - last_pos)); // add everything up to the placeholder
            size_t end = msg.find('}', i);
            if (end != std::string_view::npos) {
                std::string_view placeholder = msg.substr(i, end - i + 1);
                if (placeholder == "{iteration}") {
                    funcs.emplace_back(iteration);
                } else if (placeholder == "{iterations_total}") {
                    funcs.emplace_back(iterations_total);
                } else if (placeholder == "{chi2}") {
                    funcs.emplace_back(chi2);
                } else if (placeholder == "{chi2_best}") {
                    funcs.emplace_back(chi2_best);
                } else if (placeholder == "{chi2_penalty}") {
                    funcs.emplace_back(chi2_penalty);
                } else if (placeholder == "{chi2_no_penalty}") {
                    funcs.emplace_back(chi2_no_penalty);
                } else {
                    parts.back() += placeholder; // unknown placeholder, keep as is
                }
                i = end; // skip the placeholder
                last_pos = i+1;
            }
        }
    }
    parts.emplace_back(msg.substr(last_pos)); // add the rest of the message
    return [parts, funcs] () {
        std::string result;
        for (size_t i = 0; i < parts.size(); ++i) {
            result += parts[i];
            if (i < funcs.size()) {
                result += funcs[i]();
            }
        }
        return result;
    };
}

MessageElement::MessageElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, std::string_view message, bool log) 
    : MessageElement(owner, message, "white", log) 
{}

MessageElement::MessageElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, std::string_view message, std::string_view colour, bool log) 
    : LoopElementCallback(owner)
{
    message_func = [log, colour=std::string(colour), builder=parse_user_msg(message)] () -> void {
        if (log) {
            logging::log(builder());
        } else {
            console::print_text(builder(), console::color::parse(colour));
        }
    };
}

MessageElement::~MessageElement() = default;

void MessageElement::run() {
    assert(message_func && "MessageElement::run: message_func is not set.");
    message_func();
}