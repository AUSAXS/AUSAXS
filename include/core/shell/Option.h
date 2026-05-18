// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <string>

namespace ausaxs::shell {
    /**
     * @brief A single command-line option, consisting of a name and an associated value.
     */
    class Option {
        public:
            Option() noexcept;
            Option(const std::string& name, const std::string& value);
            virtual ~Option();

            /// @brief Get the option formatted as a string for inclusion in a command.
            std::string get() const;

            std::string name;
            std::string value;
    };

    /**
     * @brief An Option carrying a value, with convenience constructors for common value types.
     */
    class Argument : public Option {
        public:
            Argument(const std::string& name, const std::string& value);
            Argument(const std::string& name, double value);
            Argument(const std::string& name, int value);
            Argument(const std::string& name, unsigned int value);
    };

    /**
     * @brief An Option with no value, i.e. a bare command-line flag.
     */
    class Flag : public Option {
        public:
            Flag(const std::string& name);
    };
}