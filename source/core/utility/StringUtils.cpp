// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/StringUtils.h>

#include <utility/Exceptions.h>

#include <algorithm>
#include <array>
#include <cmath>

using namespace ausaxs;

std::string utility::remove_quotation_marks(std::string s) {
    if (s.size() > 1 && s[0] == '"' && s[s.size()-1] == '"') {
        return s.substr(1, s.size()-2);
    }
    return s;
}

std::string utility::remove_spaces(std::string s) {
    auto removed = std::ranges::remove(s, ' ');
    s.erase(removed.begin(), removed.end());
    return s;
}

bool utility::isdigit(char c) {
    return static_cast<bool>(std::isdigit(c));
}

bool utility::isalpha(char c) {
    return static_cast<bool>(std::isalpha(c));
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
    std::size_t i = 0;
    while (i < str.size()) {
        while (i < str.size() && str[i] == delimiter) {
            ++i;
        }
        std::size_t start = i;
        while (i < str.size() && str[i] != delimiter) {
            ++i;
        }
        if (start < i) {
            tokens.emplace_back(str.substr(start, i-start));
        }
    }
    return tokens;
}

std::vector<std::string> utility::split(std::string_view s, std::string_view delimiters) {
    std::vector<std::string> tokens;

    std::array<bool, 256> table{};
    for (auto c : delimiters) {table[c] = true;}

    // skip leading delimiters
    int start = 0;
    while (start < static_cast<int>(s.size()) && table[s[start]]) {
        ++start;
    }

    // iterate through the rest of the string
    for (int i = start; i < static_cast<int>(s.size()); ++i) {
        if (!table[s[i]]) {continue;}

        // add token to vector
        tokens.emplace_back(s.substr(start, i-start));

        ++i; // start from next char
        while (i < static_cast<int>(s.size()) && table[s[i]]) {
            ++i;
        }
        start = i;
    }

    // add last token to vector
    if (start < static_cast<int>(s.size())) {
        tokens.emplace_back(s.substr(start));
    }
    return tokens;
}

std::vector<std::string> utility::split_quoted(std::string_view s, std::string_view delimiters, char comment) {
    std::array<bool, 256> table{};
    for (auto c : delimiters) {table[c] = true;}

    std::vector<std::string> tokens;
    std::string current;
    char in_quote = 0;
    bool in_token = false;

    for (auto c : s) {
        // an unquoted comment character discards the rest of the string
        if (comment != '\0' && c == comment && in_quote == 0) {break;}

        // a quote either opens a quoted section, closes the matching one, or is a literal inside the other kind
        if (c == '"' || c == '\'') {
            if (in_quote == 0) {in_quote = c; in_token = true;}
            else if (in_quote == c) {in_quote = 0;}
            else {current += c;}
            continue;
        }

        // delimiters only delimit outside quotes
        if (in_quote == 0 && table[c]) {
            if (in_token) {
                tokens.push_back(current);
                current.clear();
                in_token = false;
            }
            continue;
        }

        current += c;
        in_token = true;
    }

    if (in_token) {tokens.push_back(current);}
    return tokens;
}

std::string utility::quote_if_needed(std::string_view s, std::string_view delimiters) {
    if (s.find_first_of(delimiters) == std::string_view::npos) {return std::string(s);}
    return '"' + std::string(s) + '"';
}

std::string utility::join(std::vector<std::string> v, std::string_view separator) {
    std::string s;
    for (int i = 0; i < static_cast<int>(v.size()); i++) {
        s += v[i];
        if (i != static_cast<int>(v.size())-1) {
            s += separator;
        }
    }
    return s;
}

std::string utility::remove_all(std::string_view s, std::string_view remove) {
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
    std::array<bool, 256> table{};
    for (auto c : remove) {table[c] = true;}

    int start = 0;
    while (start < static_cast<int>(s.size()) && table[s[start]]) {
        ++start;
    }
    return s.substr(start);
}

std::string_view utility::remove_trailing(std::string_view s, std::string_view remove) {
    std::array<bool, 256> table{};
    for (auto c : remove) {table[c] = true;}

    int end = static_cast<int>(s.size());
    while (end > 0 && table[s[end-1]]) {
        --end;
    }
    return s.substr(0, end);
}

std::string_view utility::remove_leading_and_trailing(std::string_view s, std::string_view remove) {
    std::array<bool, 256> table{};
    for (auto c : remove) {table[c] = true;}

    int start = 0;
    while (start < static_cast<int>(s.size()) && table[s[start]]) {
        ++start;
    }
    if (start == static_cast<int>(s.size())) {return s;}

    int end = static_cast<int>(s.size());
    while (end > start && table[s[end-1]]) {
        --end;
    }
    return s.substr(start, end-start);
}

std::string utility::to_lowercase(std::string_view s) {
    std::string new_s;
    for (auto c : s) {
        new_s += static_cast<char>(std::tolower(c));
    }
    return new_s;
}

bool utility::parse_bool(std::string_view s) {
    auto lower = to_lowercase(s);
    if (lower == "true" || lower == "yes" || lower == "1") {
        return true;
    }

    if (lower == "false" || lower == "no" || lower == "0") {
        return false;
    }
    throw except::invalid_argument("utility::parse_bool: \"" + std::string(s) + "\" cannot be interpreted as a boolean value.");
}

bool utility::isnumeric(std::string_view s) {
    try {std::stod(std::string(s));
    } catch (std::exception&) {return false;}
    return true;
}

bool utility::isinteger(std::string_view s) {
    try {std::stoi(std::string(s));
    } catch (std::exception&) {return false;}
    return true;
}