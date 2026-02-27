// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/ParsedArgs.h>

#include <stdexcept>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

template<>
search::ArgResult<std::string> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args, const std::string& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            return {args.at(name)[0], true};
        }
    }
    return {default_value, false};
}

template<>
search::ArgResult<std::vector<std::string>> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args, const std::vector<std::string>& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            return {args.at(name), true};
        }
    }
    return {default_value, false};
}

template<>
search::ArgResult<int> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args, const int& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            try {
                return {std::stoi(args.at(name)[0]), true};
            } catch (std::exception&) {
                throw std::invalid_argument("SequenceParser::get_arg: \"" + args.at(name)[0] + "\" cannot be interpreted as an integer.");
            }
        }
    }
    return {default_value, false};
}

template<>
search::ArgResult<std::vector<int>> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args, const std::vector<int>& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            std::vector<int> values;
            try {
                for (const auto& value : args.at(name)) {values.push_back(std::stoi(value));}
            } catch (std::exception&) {
                throw std::invalid_argument("SequenceParser::get_arg: \"" + args.at(name)[values.size()] + "\" cannot be interpreted as an integer.");
            }
            return {std::move(values), true};
        }
    }
    return {default_value, false};
}

template<>
search::ArgResult<double> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, std::vector<std::string>>& args, const double& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            try {
                return {std::stod(args.at(name)[0]), true};
            } catch (std::exception&) {
                throw std::invalid_argument("SequenceParser::get_arg: \"" + args.at(name)[0] + "\" cannot be interpreted as a decimal value.");
            }
        }
    }
    return {default_value, false};
}