#pragma once

#include <settings/SettingRef.h>

#include <string>
#include <vector>
#include <memory>

namespace settings {
    namespace io {
        /**
         * @brief By creating an instance of this class, you add another section of settings to the settings file, ready for both reading and writing.   
         * 
         * @param name The name of this section. This will be shown as a separate title line in the settings file. 
         * @param settings The settings which will be part of this section. Note that they will still be parsed if found in other sections; this is purely used for determining the write order. 
         */
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