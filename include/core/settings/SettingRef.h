// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/Type.h>
#include <utility/UtilityFwd.h>

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <stdexcept>

namespace ausaxs::settings::io::detail {
    /**
     * @brief Virtual interface for a reference to some setting. 
     */
    struct ISettingRef {
        ISettingRef(const std::vector<std::string>& name) : names(name) {}
        virtual ~ISettingRef() = default;

        /**
         * @brief Set the setting value.
         */
        virtual void set(const std::vector<std::string>&) = 0;

        /**
         * @brief Get the setting value as a string.
         */
        virtual std::string get() const = 0;

        std::vector<std::string> names; // The name of the setting.

        /**
         * @brief Get the map from setting names to value references.
         * 
         * This map contains all aliases for all settings after initialization.
         */
        static std::unordered_map<std::string, std::shared_ptr<ISettingRef>>& get_stored_settings();
    };

    /**
     * @brief A reference to a setting. 
     */
    template<typename T> struct SettingRef : public ISettingRef {
        SettingRef(T& setting, const std::vector<std::string>& names) : ISettingRef(names), settingref(setting) {}
        virtual ~SettingRef() = default;

        /**
         * @brief Set the setting value.
         */
        void set(const std::vector<std::string>&) override {
            throw std::runtime_error("settings::io::detail::SettingRef::set: missing implementation for type \"" + type(settingref) + "\".");
        }

        /**
         * @brief Get the setting value as a string.
         */
        std::string get() const override {
            throw std::runtime_error("settings::io::detail::SettingRef::get: missing implementation for type \"" + type(settingref) + "\".");
        }

        T& settingref; // A reference to the setting.
    };

    template<typename T>
    std::string type_as_string(const T&);
}

template<> std::string ausaxs::settings::io::detail::type_as_string<std::string>(const std::string&);
template<> std::string ausaxs::settings::io::detail::type_as_string<double>(const double&);
template<> std::string ausaxs::settings::io::detail::type_as_string<int>(const int&);
template<> std::string ausaxs::settings::io::detail::type_as_string<unsigned int>(const unsigned int&);
template<> std::string ausaxs::settings::io::detail::type_as_string<bool>(const bool&);
template<> std::string ausaxs::settings::io::detail::type_as_string<std::vector<std::string>>(const std::vector<std::string>&);
template<> std::string ausaxs::settings::io::detail::type_as_string<std::vector<double>>(const std::vector<double>&);
template<> std::string ausaxs::settings::io::detail::type_as_string<std::vector<int>>(const std::vector<int>&);
template<> std::string ausaxs::settings::io::detail::type_as_string<ausaxs::Limit>(const ausaxs::Limit&);

template<> std::string ausaxs::settings::io::detail::SettingRef<std::string>::get() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<double>::get() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<int>::get() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<unsigned int>::get() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<bool>::get() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<std::vector<std::string>>::get() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<std::vector<double>>::get() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<std::vector<int>>::get() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<ausaxs::Limit>::get() const;

template<> void ausaxs::settings::io::detail::SettingRef<std::string>::set(const std::vector<std::string>& str);
template<> void ausaxs::settings::io::detail::SettingRef<bool>::set(const std::vector<std::string>& str);
template<> void ausaxs::settings::io::detail::SettingRef<double>::set(const std::vector<std::string>& str);
template<> void ausaxs::settings::io::detail::SettingRef<int>::set(const std::vector<std::string>& str);
template<> void ausaxs::settings::io::detail::SettingRef<unsigned int>::set(const std::vector<std::string>& str);
template<> void ausaxs::settings::io::detail::SettingRef<std::vector<std::string>>::set(const std::vector<std::string>& str);
template<> void ausaxs::settings::io::detail::SettingRef<std::vector<double>>::set(const std::vector<std::string>& str);
template<> void ausaxs::settings::io::detail::SettingRef<std::vector<int>>::set(const std::vector<std::string>& str);
template<> void ausaxs::settings::io::detail::SettingRef<ausaxs::Limit>::set(const std::vector<std::string>& str);