// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/StringUtils.h>
#include <utility/Exceptions.h>
#include <math/ConstexprMath.h>

#include <algorithm>
#include <cmath>
#include <array>

using namespace ausaxs;

std::string utility::remove_quotation_marks(std::string s) {
    if (s.size() > 1 && s[0] == '"' && s[s.size()-1] == '"') {
        return s.substr(1, s.size()-2);
    }
    return s;
}

std::string utility::remove_spaces(std::string s) {
    std::string::iterator end_pos = std::remove(s.begin(), s.end(), ' ');
    s.erase(end_pos, s.end());
    return s;
}

std::string utility::round_double(double d, int decimals) {
    std::string s = std::to_string(std::round(d*std::pow(10, decimals))/std::pow(10, decimals));
    auto dot = s.find('.');
    if (dot == std::string::npos) {
        return s;
    }
    return s.substr(0, dot + 1 + decimals);
}

std::vector<std::string> utility::split(std::string_view str, char delimiter) {
    std::vector<std::string> tokens;
    size_t start;
    size_t i = 0;
    while (i < str.size()) {
        while (i < str.size() && str[i] == delimiter) {
            ++i;
        }
        start = i;
        while (i < str.size() && str[i] != delimiter) {
            ++i;
        }
        if (start < i) {
            tokens.push_back(std::string(str.substr(start, i-start)));
        }
    }
    return tokens;
}

std::vector<std::string> utility::split(std::string_view s, std::string_view delimiters) {
    std::vector<std::string> tokens;

    static_assert(constexpr_math::pow(2, 8*sizeof(char)) == 256, "Unexpected char size");
    std::array<bool, 256> table{};
    for (auto c : delimiters) {table[c] = true;}

    // skip leading delimiters
    unsigned int start = 0;
    while (start < s.size() && table[s[start]]) {
        ++start;
    }

    // iterate through the rest of the string
    for (unsigned int i = start; i < s.size(); ++i) {
        if (!table[s[i]]) {continue;}

        // add token to vector
        tokens.push_back(std::string(s.substr(start, i-start)));

        ++i; // start from next char
        while (i < s.size() && table[s[i]]) {
            ++i;
        }
        start = i;
    }

    // add last token to vector
    if (start < s.size()) {
        tokens.push_back(std::string(s.substr(start)));
    }
    return tokens;
}

std::string utility::join(std::vector<std::string> v, std::string_view separator) {
    std::string s;
    for (unsigned int i = 0; i < v.size(); i++) {
        s += v[i];
        if (i != v.size()-1) {
            s += separator;
        }
    }
    return s;
}

std::string utility::remove_all(std::string_view s, std::string_view remove) {
    static_assert(constexpr_math::pow(2, 8*sizeof(char)) == 256, "Unexpected char size");
    std::array<bool, 256> table{};
    for (auto c : remove) {table[c] = true;}

    std::string new_s;
    new_s.reserve(s.size());
    for (auto c : s) {
        if (!table[c]) {
            new_s.push_back(static_cast<char>(c));
        }
    }
    return new_s;
}

std::string_view utility::remove_leading(std::string_view s, std::string_view remove) {
    static_assert(constexpr_math::pow(2, 8*sizeof(char)) == 256, "Unexpected char size");
    std::array<bool, 256> table{};
    for (auto c : remove) {table[c] = true;}

    unsigned int start = 0;
    while (start < s.size() && table[s[start]]) {
        ++start;
    }
    return s.substr(start);
}

std::string_view utility::remove_trailing(std::string_view s, std::string_view remove) {
    static_assert(constexpr_math::pow(2, 8*sizeof(char)) == 256, "Unexpected char size");
    std::array<bool, 256> table{};
    for (auto c : remove) {table[c] = true;}

    unsigned int end = s.size();
    while (end > 0 && table[s[end-1]]) {
        --end;
    }
    return s.substr(0, end);
}

std::string_view utility::remove_leading_and_trailing(std::string_view s, std::string_view remove) {
    static_assert(constexpr_math::pow(2, 8*sizeof(char)) == 256, "Unexpected char size");
    std::array<bool, 256> table{};
    for (auto c : remove) {table[c] = true;}

    unsigned int start = 0;
    while (start < s.size() && table[s[start]]) {
        ++start;
    }
    if (start == s.size()) {return s;}

    unsigned int end = s.size();
    while (end > start && table[s[end-1]]) {
        --end;
    }
    return s.substr(start, end-start);
}

std::string utility::to_lowercase(std::string_view s) {
    std::string new_s;
    for (auto c : s) {
        new_s += std::tolower(c);
    }
    return new_s;
}

bool utility::parse_bool(std::string_view s) {
    auto lower = to_lowercase(s);
    if (lower == "true" || lower == "yes" || lower == "1") {
        return true;
    } else if (lower == "false" || lower == "no" || lower == "0") {
        return false;
    }
    throw except::invalid_argument("utility::parse_bool: \"" + std::string(s) + "\" cannot be interpreted as a boolean value.");
}
