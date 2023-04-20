#pragma once

#include <string>
#include <vector>
#include <unordered_map>

namespace settings {
    namespace detail {
        /**
         * @brief A smart options class. This is used to store settings in a way that allows them to be read from and written to a file.
         */
        struct ISmartOption {
            virtual ~ISmartOption() = default;

            /**
             * @brief Set this option to the given value.
             */
            virtual void set(const std::vector<std::string>& value) = 0;

            /**
             * @brief Get the current value of this option as a std::string.
             */
            virtual std::string get() const = 0;

            /**
             * @brief A map of all options, indexed by their name.
             * @note This is used to find the correct option when parsing a settings file.
             */
            inline static std::unordered_map<std::string, ISmartOption*> all_options;
        };

        /**
         * @brief A smart option for a specific type.
         * @tparam T The type of the option.
         */
        template<typename T>
        struct SmartOption : public ISmartOption {
            /**
             * @brief Construct a new option.
             * 
             * @param value The initial value of the option.
             * @param name Accepted name for the option. These will be written to and can be read from files.
             */
            SmartOption(const T& value, const std::initializer_list<std::string>& name) : value(value) {
                for (const auto& n : name) {all_options.emplace(n, this);}
            }

            /**
             * @brief Construct a new option.
             * 
             * @param value The initial value of the option.
             * @param name Accepted name for the option. This will be written to and can be read from files.
             */
            SmartOption(const T& value, const std::string& name) : value(value) {
                all_options.emplace(name, this);
            }
            
            /**
             * @brief Construct a new unnamed option.
             * @note This option will not be accessible by name and cannot interact with files.
             * 
             * @param value The initial value of the option.
             */
            SmartOption(const T& value) : value(value) {}

            ~SmartOption() override = default;

            /**
             * @brief Get the current value of this option.
             */
            [[nodiscard]] operator T() const {
                return value;
            }

            T& operator=(const T& other) {
                value = other;
                return value;
            }

            friend T operator+(const SmartOption& lhs, const T& rhs) {return lhs.value + rhs;}

            /**
             * @brief Print this option to the given stream.
             */
            friend std::ostream& operator<<(std::ostream& os, const SmartOption& option) {
                os << option.value;
                return os;
            }

            /**
             * @brief Set this option to the given value.
             */
            void set(const std::vector<std::string>& value) override;

            /**
             * @brief Get the current value of this option as a string.
             */
            [[nodiscard]] std::string get() const override;

            T value;
        };
    }
}