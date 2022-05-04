#pragma once

#include <string>
#include <sstream>

namespace utility {
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
     *        Note that the input string is modified. 
     */
    std::string remove_spaces(std::string s);

    template<typename T>
    T extract_number(std::string s) {
        std::stringstream ss(extract_number<std::string>(s));
        T val; ss >> val;
        return val;
    }

    /**
     * @brief Remove the extension from a filename. 
     *        This is just a simple wrapper around filesystem::path::replace_extension.
     */
    std::string remove_extension(std::string path);

    /**
     * @brief Create all parent directories of the path.
     */
    void create_directories(std::string& path);

    /**
     * @brief Print a warning message. The text will be red in the terminal. 
     */
    void print_warning(std::string text);
}