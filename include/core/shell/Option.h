// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <string>

namespace ausaxs::shell {
    class Option {
        public:
            Option() noexcept;
            Option(const std::string& name, const std::string& value);
            virtual ~Option();

            std::string get() const;

            std::string name;
            std::string value;            
    };

    class Argument : public Option {
        public:
            Argument(const std::string& name, const std::string& value);
            Argument(const std::string& name, double value);
            Argument(const std::string& name, int value);
            Argument(const std::string& name, unsigned int value);
    };

    class Flag : public Option {
        public:
            Flag(const std::string& name);
    };
}