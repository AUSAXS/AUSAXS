// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <string>
#include <vector>
#include <unordered_map>

namespace ausaxs::rigidbody::sequencer::detail {
    template<class Args>
    std::vector<std::string> get_arg_names(const std::unordered_map<Args, std::vector<std::string>>& args_map) {
        std::vector<std::string> args; 
        args.reserve(args_map.size() * 2);
        for (const auto& [_, arg_names] : args_map) {
            for (const auto& name : arg_names) {
                args.emplace_back(name);
            }
        }
        return args;
    }
}