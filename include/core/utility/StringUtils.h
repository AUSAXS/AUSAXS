#pragma once

#include <string_view>
#include <string>
#include <vector>

namespace ausaxs::utility {
    /**
     * @brief Remove spaces from both ends of a string. 
     */
    std::string remove_spaces(std::string s);

    /**
     * @brief Remove quotation marks from both ends of a string. 
     *        If the string is not enclosed in quotation marks, it is returned unchanged.
     */
    std::string remove_quotation_marks(std::string s);

    /**
     * @brief Convert a string to lowercase.
     */
    std::string to_lowercase(std::string_view s);

    /**
     * @brief Round a double to a given number of decimal places.
     */
    std::string round_double(double d, int decimals);

    /**
     * @brief Split a string at a given delimiter.
     *        Consecutive delimiters are treated as a single delimiter. 
     */
    std::vector<std::string> split(std::string_view s, char delimiter);

    /**
     * @brief Split a string at the given delimiters.
     *        Consecutive delimiters are treated as a single delimiter. 
     */
    std::vector<std::string> split(std::string_view s, std::string_view delimiters);

    /**
     * @brief Join a vector of strings into a single string. The separator will be inserted after each element except the last. 
     */
    std::string join(std::vector<std::string> v, std::string_view separator);

    /**
     * @brief Remove all occurrences of the characters in 'remove' from the string. 
     */
    std::string remove_all(std::string_view s, std::string_view remove);

    /**
     * @brief Parse a string as a boolean value.
     *        The following strings are considered true: "true", "yes", "1".
     *        The following strings are considered false: "false", "no", "0".
     *        Other input will throw an exception.
     */
    bool parse_bool(std::string_view s);
}