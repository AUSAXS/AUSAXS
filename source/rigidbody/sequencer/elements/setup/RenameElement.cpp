// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/RenameElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <utility/observer_ptr.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

RenameElement::RenameElement(observer_ptr<Sequencer> owner, std::string_view old_name, std::string_view new_name) {
    auto& body_names = owner->setup()._get_body_names();
    auto it = body_names.find(std::string{old_name});
    unsigned int index = it->second;
    body_names.erase(it);
    body_names.emplace(std::string{new_name}, index);
}

RenameElement::~RenameElement() = default;

void RenameElement::run() {}

std::vector<std::string> RenameElement::_valid_arguments() {
    return {};
}

std::unique_ptr<GenericElement> RenameElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("rename", "Unexpected named argument \"" + args.named.begin()->first + "\".");}
    if (args.inlined.size() != 2) {throw except::parse_error(
        "rename", "Invalid number of inline arguments. Expected [old name] [new name], but got " + std::to_string(args.inlined.size()) + "."
    );}

    const auto& body_names = owner->_get_sequencer()->setup()._get_body_names();
    std::string old_name = args.inlined[0];
    std::string new_name = args.inlined[1];
    if (!body_names.contains(old_name)) {throw except::parse_error("rename", "Body name \"" + old_name + "\" not found.");}
    if (old_name == new_name) {throw except::parse_error("rename", "Cannot rename body \"" + old_name + "\" to itself.");}
    if (body_names.contains(new_name)) {throw except::parse_error("rename", "Body name \"" + new_name + "\" already exists.");}

    return std::make_unique<RenameElement>(owner->_get_sequencer(), old_name, new_name);
}
