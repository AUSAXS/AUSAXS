// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/SettingsIORegistry.h>
#include <utility/Exceptions.h>
#include <utility/observer_ptr.h>

using namespace ausaxs;
using namespace ausaxs::settings::io;

std::vector<observer_ptr<SettingSection>>& SettingSection::get_sections() {
    static std::vector<observer_ptr<SettingSection>> sections;
    return sections;
}

SettingSection::SettingSection(std::string_view name, std::initializer_list<std::shared_ptr<detail::ISettingRef>> settings) : name(name), settings(settings) {
    auto& stored_settings = detail::ISettingRef::get_stored_settings();
    for (auto& setting : settings) {
        for (auto& name : setting->names) {
            if (stored_settings.contains(name)) {
                throw std::runtime_error("Settings::add: Duplicate setting name: \"" + name + "\".");
            }
            stored_settings[name] = std::move(setting);
        }
    }
    get_sections().push_back(this);
}