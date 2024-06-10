/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/utility/Utility.h>
#include <md/utility/ConsoleColor.h>

#include <sstream>

std::vector<std::string> utility::split(std::string str, char delimiter) {
    std::string token;
    std::stringstream ss(str);
    std::vector<std::string> tokens;
    while(std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::string> utility::split(std::string str, std::string delimiters) {
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

std::string utility::join(const std::vector<std::string>& tokens, const std::string& delimiter) {
    std::string str;
    for (unsigned int i = 0; i < tokens.size(); i++) {
        str += tokens[i];
        if (i != tokens.size() - 1) {
            str += delimiter;
        }
    }
    return str;
}

void utility::print_warning(const std::string& text) {
    console::print(text, console::color::red);
}

void utility::print_success(const std::string& text) {
    console::print(text, console::color::green);
}

void utility::print_info(const std::string& text) {
    console::print(text, console::color::lightblue);
}
