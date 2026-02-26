// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <vector>
#include <string>
#include <unordered_map>

namespace ausaxs::rigidbody::sequencer {
    struct ParsedArgs {
        using Key = std::string;
        using Value = std::string;
        struct Args {
            struct Arg {
                int line_number;
                Value value;
            };

            int line_number;
            std::vector<Arg> args;
        };

        struct InlineArgs {
            int line_number;
            std::vector<Value> values;
        };

        InlineArgs inlined;
        std::unordered_map<Key, Args> named;
    };
}