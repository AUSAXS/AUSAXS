/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <utility/Utility.h>
#include <utility/Console.h>

using namespace ausaxs;

bool utility::approx(double v1, double v2, double abs, double eps) {
    return std::abs(v1 - v2) <= std::max(abs, eps * std::max(std::abs(v1), std::abs(v2)));
}

bool utility::equal(double a, double b, double c) {
    return a == b && b == c;
}

std::string utility::uid() {
    static unsigned int i = 0;
    return std::to_string(i++);
}

std::string utility::uid(const std::string& s) {return s + uid();}

std::ostream& utility::detail::operator<<(std::ostream& os, const __dummy& obj) {
    os << obj.s;
    return os;
}

std::string utility::round(double val, unsigned int decimals) noexcept {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(decimals) << val;
    return ss.str();
}

utility::detail::__dummy utility::fixedwidth(double number, unsigned int width) {
    std::string s = std::to_string(number);
    // remove unnecessary 0s
    unsigned int end = s.size();
    for (unsigned int i = end; i > 0; i--) {
        if (s[i-1] != '0') {
            end = i+1;
            break;
        }
    }
    if (end < s.size()) {
        s.resize(end);
    }

    std::string o;
    for (unsigned int i = 0; i < width; i++) {
        if (i < s.size()) {
            o += s[i];
        } else {
            o += ' ';
        }
    }

    // check how lossy the conversion was
    #ifdef DEBUG
        double d = std::stod(o);
        if (!approx(d, number, 1e-3)) {
            console::print_warning("Fixed-width conversion of " + std::to_string(number) + " to " + o + " is lossy.");
        }
    #endif
    
    return {std::move(o)};
} 