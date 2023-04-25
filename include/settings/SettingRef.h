#pragma once

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <utility/Type.h>

namespace settings {
    namespace io {
        namespace detail {
            /**
             * @brief Virtual interface for a reference to some setting. 
             */
            struct ISettingRef {
                ISettingRef(const std::vector<std::string>& name) : names(name) {}
                virtual ~ISettingRef() = default;

                /**
                 * @brief Set the setting value.
                 */
                virtual void set(const std::vector<std::string>&) {
                    throw std::runtime_error("settings::io::detail::ISettingRef::set: This should never happen.");
                }

                /**
                 * @brief Get the setting value as a string.
                 */
                virtual std::string get() const {
                    throw std::runtime_error("settings::io::detail::ISettingRef::get: This should never happen.");
                }

                std::vector<std::string> names; // The name of the setting.
                inline static std::unordered_map<std::string, std::shared_ptr<ISettingRef>> stored_settings;
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

        }
    }
}