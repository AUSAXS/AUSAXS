// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/SettingsElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <settings/SettingRef.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

SettingsElement::SettingsElement(observer_ptr<Sequencer> owner, std::string name, std::string value)
    : owner(owner), name(std::move(name)), value(std::move(value)) 
{
    auto& settings = settings::io::detail::ISettingRef::get_stored_settings();
    settings.at(this->name)->set({this->value});
}

SettingsElement::~SettingsElement() = default;

void SettingsElement::run() {}

std::vector<std::string> SettingsElement::_valid_arguments() {
    return {};
}

InlineSignature SettingsElement::_valid_inline_arguments() {
    return {.names = {"setting", "value"}, .min = 2, .max = 2};
}

std::unique_ptr<GenericElement> SettingsElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    const std::string name = args.inlined[0];
    const auto& settings = settings::io::detail::ISettingRef::get_stored_settings();
    if (!settings.contains(name)) {
        throw except::parse_error("settings", "Unknown setting \"" + name + "\".");
    }

    return std::make_unique<SettingsElement>(owner->_get_sequencer(), name, args.inlined[1]);
}
