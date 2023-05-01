#pragma once

#include <string>
#include <vector>

namespace utility {
    /**
     * @brief Remove spaces from both ends of a string. 
     */
    std::string remove_spaces(std::string s);

    /**
     * @brief Remove quotation marks from both ends of a string. 
     */
    std::string remove_quotation_marks(std::string s);

    /**
     * @brief Convert a string to lowercase.
     */
    std::string to_lowercase(const std::string& s);

    /**
     * @brief Split a string at a given delimiter.
     *        Consecutive delimiters are treated as a single delimiter. 
     */
    std::vector<std::string> split(const std::string& s, char delimiter);

    /**
     * @brief Split a string at the given delimiters.
     *        Consecutive delimiters are treated as a single delimiter. 
     */
    std::vector<std::string> split(const std::string& s, const std::string& delimiters);

    /**
     * @brief Join a vector of strings into a single string. The separator will be inserted after each element except the last. 
     */
    std::string join(std::vector<std::string> v, const std::string& separator);

    /**
     * @brief Remove all occurrences of the characters in 'remove' from the string. 
     */
    std::string remove_all(const std::string& s, const std::string& remove);
}