#pragma once

#include <string>

namespace shell {
    class Option {
        public:
            Option() = default;
            Option(const std::string& name, const std::string& value) : name(name), value(value) {}
            virtual ~Option() = default;

            std::string get() const {
                if (value.empty()) {
                    return name;
                }
                return name + " " + value;
            }

            std::string name;
            std::string value;            
    };

    class Argument : public Option {
        public:
            Argument(const std::string& name, const std::string& value) : Option(name, value) {}
            Argument(const std::string& name, double value) : Option(name, std::to_string(value)) {}
            Argument(const std::string& name, int value) : Option(name, std::to_string(value)) {}
            Argument(const std::string& name, unsigned int value) : Option(name, std::to_string(value)) {}
    };

    class Flag : public Option {
        public:
            Flag(const std::string& name) : Option(name, "") {}
    };
}