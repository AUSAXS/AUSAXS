#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <stdexcept>

#include <iostream>
namespace md {
    struct ISmartOption {
        virtual ~ISmartOption() = default;
        virtual void set(const std::string& value) = 0;
        virtual std::string get() const = 0;
        static std::unordered_map<std::string, ISmartOption*>& get_all_options();
    };

    template<typename T>
    struct SmartOption : public ISmartOption {
        SmartOption(const T& value, const std::vector<std::string>& name) : value(value) {
            for (const auto& n : name) {get_all_options().emplace(n, this);}
        }
        ~SmartOption() override = default;

        operator T() const {
            return value;
        }

        friend std::ostream& operator<<(std::ostream& os, const SmartOption& option) {
            os << option.value;
            return os;
        }

        void set(const std::string& value) override;
        std::string get() const override;

        T value;
    };

    struct settings {
        static bool is_comment_char(char c);
        static void parse_option(const std::string& name, const std::string& value);
        static void read(const std::string& filename);
        static void write(const std::string& filename);
        static bool discover(std::string path);
    };
}