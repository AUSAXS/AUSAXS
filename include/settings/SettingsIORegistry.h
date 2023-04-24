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

// template<> std::string settings::io::detail::SettingRef<std::string>::get() const;
// template<> std::string settings::io::detail::SettingRef<double>::get() const;
// template<> std::string settings::io::detail::SettingRef<int>::get() const;
// template<> std::string settings::io::detail::SettingRef<unsigned int>::get() const;
// template<> std::string settings::io::detail::SettingRef<bool>::get() const;
// template<> std::string settings::io::detail::SettingRef<std::vector<std::string>>::get() const;
// template<> std::string settings::io::detail::SettingRef<std::vector<double>>::get() const;
// template<> std::string settings::io::detail::SettingRef<std::vector<int>>::get() const;

// template<> void settings::io::detail::SettingRef<std::string>::set(const std::vector<std::string>& str);
// template<> void settings::io::detail::SettingRef<bool>::set(const std::vector<std::string>& str);
// template<> void settings::io::detail::SettingRef<double>::set(const std::vector<std::string>& str);
// template<> void settings::io::detail::SettingRef<int>::set(const std::vector<std::string>& str);
// template<> void settings::io::detail::SettingRef<unsigned int>::set(const std::vector<std::string>& str);
// template<> void settings::io::detail::SettingRef<std::vector<std::string>>::set(const std::vector<std::string>& str);
// template<> void settings::io::detail::SettingRef<std::vector<double>>::set(const std::vector<std::string>& str);
// template<> void settings::io::detail::SettingRef<std::vector<int>>::set(const std::vector<std::string>& str);