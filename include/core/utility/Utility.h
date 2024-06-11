#pragma once

#include <string>
#include <sstream>
#include <iomanip>

namespace utility {
    /**
     * @brief Print an element with the given width.
     */
    template<typename T>
    struct print_element {
        print_element(const T& t, int width) : t(t), width(width) {}

        friend std::ostream& operator<<(std::ostream& os, const print_element<T> e) noexcept {
            std::stringstream ss; ss << e.t;
            std::string val = ss.str();
            if (val.size() > e.width) {val = val.substr(0, e.width);}

            os << std::left << std::setw(e.width) << e.t; return os;
        }

        T t;
        unsigned int width;
    };

    /**
     * @brief Check if two numbers are approximately equal. 
     * 
     * @param v1 First value.
     * @param v2 Second value. 
     * @param abs Absolute tolerance. 
     * @param eps Relative tolerance. 
     */
    bool approx(double v1, double v2, double abs = 1e-6, double eps = 0.01);

    /**
     * @brief Check if three values are equal.
     */
    bool equal(double a, double b, double c);

    /**
     * @brief Get a unique identifier.
     */
    std::string uid();

    /**
     * @brief Append a unique identifier to a string.
     */
    std::string uid(const std::string& s);

    /**
     * @brief Round a number to a string.
     */
    std::string round(double val, unsigned int decimals) noexcept;

    namespace detail {
        // Dummy object for fixed-length printing of numbers. 
        // std::setprecision does *not* count leading zeros, which breaks our strict formatting.
        struct __dummy {
            std::string s;
        };

        std::ostream& operator<<(std::ostream& os, const __dummy& obj);
    }

    detail::__dummy fixedwidth(double number, unsigned int width);
}