// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/SettingsHelper.h>
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

        /**
         * @brief Get the type of the setting as a string. This is primarily used for introspection in the Python wrapper. 
         */
        virtual std::string type() const = 0;

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

        void set(const std::vector<std::string>&) override {
            throw std::runtime_error("settings::io::detail::SettingRef::set: missing implementation for type \"" + ausaxs::type(settingref) + "\".");
        }

        std::string get() const override {
            throw std::runtime_error("settings::io::detail::SettingRef::get: missing implementation for type \"" + ausaxs::type(settingref) + "\".");
        }

        std::string type() const override {
            throw std::runtime_error("settings::io::detail::SettingRef::type: missing implementation for type \"" + ausaxs::type(settingref) + "\".");
        }

        T& settingref; // A reference to the setting.
    };

    // Specialization for Setting<T> wrapper
    template<typename T>
    struct SettingRef<ausaxs::settings::detail::Setting<T>> : public ISettingRef {
        SettingRef(ausaxs::settings::detail::Setting<T>& setting, const std::vector<std::string>& names)
            : ISettingRef(names), settingref(setting) {}
        virtual ~SettingRef() = default;

        void set(const std::vector<std::string>& values) override {
            // Delegate to SettingRef<T>
            SettingRef<T> ref(settingref.value, names);
            ref.set(values);
        }

        std::string get() const override {
            // Delegate to SettingRef<T>
            SettingRef<T> ref(const_cast<T&>(settingref.value), names);
            return ref.get();
        }

        std::string type() const override {
            // Delegate to SettingRef<T>
            SettingRef<T> ref(const_cast<T&>(settingref.value), names);
            return ref.type();
        }

        ausaxs::settings::detail::Setting<T>& settingref;
    };
}

template<> std::string ausaxs::settings::io::detail::SettingRef<std::string>::type() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<double>::type() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<int>::type() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<unsigned int>::type() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<bool>::type() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<std::vector<std::string>>::type() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<std::vector<double>>::type() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<std::vector<int>>::type() const;
template<> std::string ausaxs::settings::io::detail::SettingRef<ausaxs::Limit>::type() const;

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