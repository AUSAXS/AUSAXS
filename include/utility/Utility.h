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
     * @brief Append a string to the stem of a path. 
     *        Example: path = "dir/file.txt", s = "_raw" --> "dir/file_raw.txt"
     */
    std::string stem_append(std::string path, std::string s);

    /**
     * @brief Get the stem of a path.
     */
    std::string stem(std::string path);

    /**
     * @brief Check if three values are equal.
     */
    bool equal(double a, double b, double c);

    /**
     * @brief Print a warning message. The text will be red in the terminal. 
     */
    void print_warning(std::string text);

    /**
     * @brief Get a unique identifier.
     */
    std::string uid();

    /**
     * @brief Append a unique identifier to a string.
     */
    std::string uid(std::string s);
}