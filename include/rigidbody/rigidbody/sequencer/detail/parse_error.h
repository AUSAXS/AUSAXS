#pragma once

#include <rigidbody/sequencer/detail/ParsedArgs.h>

#include <stdexcept>
#include <string>
#include <numeric>

namespace ausaxs::rigidbody::sequencer::except {
    struct parse_error : public std::invalid_argument {
        explicit parse_error(std::string_view element, std::string_view msg) : std::invalid_argument(
            "Error parsing element \"" + std::string(element) + "\": " + std::string(msg)
        ) {}

        explicit parse_error(std::string_view element, const ParsedArgs::Args::Arg& arg, std::string_view msg) : std::invalid_argument(
            "Error parsing element \"" + std::string(element) + "\" at line " + std::to_string(arg.line_number) + ": "
            "Invalid argument \"" + arg.value + "\".\n\t" + std::string(msg)
        ) {}

        explicit parse_error(std::string_view element, const ParsedArgs::InlineArgs& args, std::string_view msg) : std::invalid_argument(
            "Error parsing element \"" + std::string(element) + "\" at line " + std::to_string(args.line_number) + "."
            "Invalid inline arguments: " + std::accumulate(args.values.begin(), args.values.end(), std::string(),
                [] (const std::string& a, const std::string& b) {
                    return a.empty() ? b : a + ", " + b;
                }
            ) + ".\n\t" + std::string(msg)
        ) {}
    };
}