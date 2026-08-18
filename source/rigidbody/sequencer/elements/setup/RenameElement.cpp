// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/RenameElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <utility/observer_ptr.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

RenameElement::RenameElement(observer_ptr<Sequencer> owner, std::string_view old_name, std::string_view new_name) {
    owner->setup()._body_name_registry().rename(old_name, new_name);
}

RenameElement::~RenameElement() = default;

void RenameElement::run() {}

std::vector<std::string> RenameElement::_valid_arguments() {
    return {};
}

InlineSignature RenameElement::_valid_inline_arguments() {
    return {.names = {"old name", "new name"}, .min = 2, .max = 2};
}

std::unique_ptr<GenericElement> RenameElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    const auto& body_names = owner->_get_sequencer()->setup()._body_name_registry();
    std::string old_name = args.inlined[0];
    std::string new_name = args.inlined[1];
    if (!body_names.contains(old_name)) {throw except::parse_error("rename", "Body name \"" + old_name + "\" not found.");}
    if (old_name == new_name) {throw except::parse_error("rename", "Cannot rename body \"" + old_name + "\" to itself.");}
    if (body_names.contains(new_name)) {throw except::parse_error("rename", "Body name \"" + new_name + "\" already exists.");}

    return std::make_unique<RenameElement>(owner->_get_sequencer(), old_name, new_name);
}
