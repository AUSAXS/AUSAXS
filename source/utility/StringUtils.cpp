#include <utility/StringUtils.h>

#include <sstream>

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

std::vector<std::string> utility::split(const std::string& str, char delimiter) {
    std::string token;
    std::stringstream ss(str);
    std::vector<std::string> tokens;
    while(std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::string> utility::split(const std::string& str, const std::string& delimiters) {
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
        tokens.push_back(str.substr(start, i-start));
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
        tokens.push_back(str.substr(start));
    }
    return tokens;
}

std::string utility::join(std::vector<std::string> v, const std::string& separator) {
    std::string s;
    for (unsigned int i = 0; i < v.size(); i++) {
        s += v[i];
        if (i != v.size()-1) {
            s += separator;
        }
    }
    return s;
}

std::string utility::remove_all(const std::string& s, const std::string& remove) {
    std::string new_s;
    for (auto c : s) {
        if (remove.find(c) == std::string::npos) {
            new_s += c;
        }
    }
    return new_s;
}

std::string utility::to_lowercase(const std::string& s) {
    std::string new_s;
    for (auto c : s) {
        new_s += std::tolower(c);
    }
    return new_s;
}