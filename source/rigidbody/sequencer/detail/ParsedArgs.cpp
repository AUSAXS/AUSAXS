// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/ParsedArgs.h>

#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

template<>
search::ArgResult<std::string> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, ParsedArgs::Args>& args, const std::string& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            return {.value=args.at(name)[0].str, .found=true};
        }
    }
    return {.value=default_value, .found=false};
}

template<>
search::ArgResult<std::vector<std::string>> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, ParsedArgs::Args>& args, const std::vector<std::string>& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            std::vector<std::string> result;
            for (const auto& arg : args.at(name).args) {result.push_back(arg.str);}
            return {.value=std::move(result), .found=true};
        }
    }
    return {.value=default_value, .found=false};
}

template<>
search::ArgResult<int> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, ParsedArgs::Args>& args, const int& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            try {
                return {.value=std::stoi(args.at(name)[0].str), .found=true};
            } catch (std::exception&) {
                throw except::invalid_argument("SequenceParser::get_arg: \"" + args.at(name)[0].str + "\" cannot be interpreted as an integer.");
            }
        }
    }
    return {.value=default_value, .found=false};
}

template<>
search::ArgResult<std::vector<int>> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, ParsedArgs::Args>& args, const std::vector<int>& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            std::vector<int> values;
            try {
                for (const auto& arg : args.at(name).args) {values.push_back(std::stoi(arg.str));}
            } catch (std::exception&) {
                throw except::invalid_argument("SequenceParser::get_arg: \"" + args.at(name).args[values.size()].str + "\" cannot be interpreted as an integer.");
            }
            return {.value=std::move(values), .found=true};
        }
    }
    return {.value=default_value, .found=false};
}

template<>
search::ArgResult<double> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, ParsedArgs::Args>& args, const double& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            try {
                return {.value=std::stod(args.at(name)[0].str), .found=true};
            } catch (std::exception&) {
                throw except::invalid_argument("SequenceParser::get_arg: \"" + args.at(name)[0].str + "\" cannot be interpreted as a decimal value.");
            }
        }
    }
    return {.value=default_value, .found=false};
}