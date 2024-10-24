/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <utility/StringUtils.h>

#include <algorithm>
#include <sstream>

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

std::vector<std::string> utility::split(std::string_view str, char delimiter) {
    std::string token;
    std::stringstream ss((std::string(str)));
    std::vector<std::string> tokens;
    while(std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::string> utility::split(std::string_view str, std::string_view delimiters) {
    std::vector<std::string> tokens;

    auto is_delimiter = [&delimiters] (char c) {
        return delimiters.find(c) != std::string::npos;
    };

    // skip leading delimiters
    unsigned int start = 0;
    for (; start < str.size(); start++) {
        if (!is_delimiter(str[start])) {
            break;
        }
    }

    // iterate through the rest of the string
    for (unsigned int i = start; i < str.size(); i++) {
        if (!is_delimiter(str[i])) {
            continue;
        }

        // add token to vector
        tokens.push_back(std::string(str.substr(start, i-start)));
        start = ++i;

        // skip consecutive delimiters
        for (; start < str.size(); start++) {
            if (!is_delimiter(str[start])) {
                break;
            }
        }
        i = start;
    }

    // add last token to vector
    if (start < str.size()) {
        tokens.push_back(std::string(str.substr(start)));
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
    std::string new_s;
    for (auto c : s) {
        if (remove.find(c) == std::string::npos) {
            new_s += c;
        }
    }
    return new_s;
}

std::string utility::to_lowercase(std::string_view s) {
    std::string new_s;
    for (auto c : s) {
        new_s += std::tolower(c);
    }
    return new_s;
}