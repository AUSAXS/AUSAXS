#include <settings/SettingsIORegistry.h>
#include <utility/Exceptions.h>

using namespace settings::io;

std::vector<SettingSection> SettingSection::sections = std::vector<SettingSection>();

SettingSection::SettingSection(std::string_view name, std::initializer_list<std::shared_ptr<detail::ISettingRef>> settings) : name(name), settings(settings) {
    for (auto& setting : settings) {
        for (auto& name : setting->names) {
            if (detail::ISettingRef::stored_settings.contains(name)) {
                throw std::runtime_error("Settings::add: Duplicate setting name: \"" + name + "\".");
            }
            detail::ISettingRef::stored_settings[name] = std::move(setting);
        }
    }
}