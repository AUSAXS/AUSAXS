#pragma once

#include <utility/settings/SettingRef.h>

#include <string>
#include <vector>
#include <memory>

namespace settings {
    namespace io {
        struct SettingSection {
            SettingSection(std::string_view name, std::initializer_list<std::shared_ptr<detail::ISettingRef>> settings);

            std::string name;
            std::vector<std::shared_ptr<detail::ISettingRef>> settings;
        };

        template<typename T> std::unique_ptr<detail::ISettingRef> create(T& setting, const std::string& name) {
            detail::SettingRef temp(setting, {name});
            return std::make_unique<detail::ISettingRef>(std::move(temp));
            // return std::make_unique<detail::SettingRef<T>>(setting, {name});
        }

        template<typename T> std::unique_ptr<detail::ISettingRef> create(T& setting, const std::initializer_list<std::string>& names) {
            detail::SettingRef temp(setting, names);
            return std::make_unique<detail::ISettingRef>(std::move(temp));
            // return std::make_unique<detail::SettingRef<T>>(setting, names);
        }

        namespace detail {
            extern std::vector<SettingSection> settings_sections;
        }
    }
}