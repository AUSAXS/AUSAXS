#pragma once

#include <settings/SettingRef.h>

#include <string>
#include <vector>
#include <memory>

namespace settings {
    namespace io {
        struct SettingSection {
            SettingSection(std::string_view name, std::initializer_list<std::shared_ptr<detail::ISettingRef>> settings);

            std::string name;
            std::vector<std::shared_ptr<detail::ISettingRef>> settings;
            inline static std::vector<SettingSection> sections;
        };

        template<typename T> std::unique_ptr<detail::SettingRef<T>> create(T& setting, const std::string& name) {
            return std::make_unique<detail::SettingRef<T>>(setting, std::vector<std::string>({name}));
        }

        template<typename T> std::unique_ptr<detail::SettingRef<T>> create(T& setting, const std::initializer_list<std::string>& names) {
            return std::make_unique<detail::SettingRef<T>>(setting, names);
        }
    }
}