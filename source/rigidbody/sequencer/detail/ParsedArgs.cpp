// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/ParsedArgs.h>

#include <utility/Exceptions.h>
#include <utility/StringUtils.h>

#include <charconv>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

int search::parse_int(std::string_view token) {
    std::string_view s = token.starts_with('+') ? token.substr(1) : token;
    int value;
    auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
    if (s.empty() || ec != std::errc() || end != s.data() + s.size()) {
        throw except::invalid_argument("SequenceParser: \"" + std::string(token) + "\" cannot be interpreted as an integer.");
    }
    return value;
}

double search::parse_double(std::string_view token) {
    std::string_view s = token.starts_with('+') ? token.substr(1) : token;
    std::string_view digits = s.starts_with('-') ? s.substr(1) : s;
    bool spelled_finite = !digits.empty() && (utility::isdigit(digits[0]) || digits[0] == '.');

    double value;
    auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
    if (!spelled_finite || ec != std::errc() || end != s.data() + s.size()) {
        throw except::invalid_argument("SequenceParser: \"" + std::string(token) + "\" cannot be interpreted as a finite decimal value.");
    }
    return value;
}

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
            return {.value=parse_int(args.at(name)[0].str), .found=true};
        }
    }
    return {.value=default_value, .found=false};
}

template<>
search::ArgResult<std::vector<int>> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, ParsedArgs::Args>& args, const std::vector<int>& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            std::vector<int> values;
            for (const auto& arg : args.at(name).args) {values.push_back(parse_int(arg.str));}
            return {.value=std::move(values), .found=true};
        }
    }
    return {.value=default_value, .found=false};
}

template<>
search::ArgResult<double> search::get_arg(std::vector<std::string>& names, const std::unordered_map<std::string, ParsedArgs::Args>& args, const double& default_value) {
    for (const auto& name : names) {
        if (args.contains(name)) {
            return {.value=parse_double(args.at(name)[0].str), .found=true};
        }
    }
    return {.value=default_value, .found=false};
}