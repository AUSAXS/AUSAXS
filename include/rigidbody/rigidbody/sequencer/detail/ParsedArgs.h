// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <vector>
#include <string>
#include <unordered_map>

namespace ausaxs::rigidbody::sequencer {
    namespace search {
        template<typename T>
        struct ArgResult {T value; bool found;};

        template<typename T>
        ArgResult<T> get_arg(std::vector<std::string>& valid_keys, const std::unordered_map<std::string, std::vector<std::string>>& args, const T& default_value = T());
    }

    struct ParsedArgs {
        using Key = std::string;
        using Value = std::string;
        struct Args {
            std::size_t size() const {return args.size();}
            bool empty() const {return args.empty();}

            struct Arg {
                int line_number;
                Value value;
            };

            int line_number;
            std::vector<Arg> args;
        };

        struct InlineArgs {
            std::size_t size() const {return values.size();}
            bool empty() const {return values.empty();}

            int line_number;
            std::vector<Value> values;
        };

        template<typename T>
        search::ArgResult<T> get(std::vector<Key>& valid_keys, const T& default_value = T()) const;

        InlineArgs inlined;
        std::unordered_map<Key, Args> named;
    };
}

template<typename T>
ausaxs::rigidbody::sequencer::search::ArgResult<T> ausaxs::rigidbody::sequencer::ParsedArgs::get(
    std::vector<std::string>& valid_keys, const T& default_value
) const {
    return search::get_arg(valid_keys, named, default_value);
}