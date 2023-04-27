#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

//! Split some of this functionality into the File class.
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