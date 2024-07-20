#pragma once

#include <string_view>

namespace console {
    /**
     * @brief Print a warning message. The text will be red in the terminal. 
     */
    void print_warning(std::string_view text);

    /**
     * @brief Print a success message. The text will be green in the terminal. 
     */
    void print_success(std::string_view text);

    /**
     * @brief Print a info message. The text will be blue in the terminal. 
     *        Should only be used as a header for a info section. Use tabs to indent other text in the section. 
     */
    void print_info(std::string_view text);

    /**
     * @brief Print a text message. The text will be white in the terminal. 
     */
    void print_text(std::string_view text);

    /**
     * @brief Print a minor text message. The text will be white in the terminal. 
     *        Will only be printed if minor messages are enabled.
     */
    void print_text_minor(std::string_view text);

    /**
     * @brief Add another indentation to print_text output.
     */
    void indent();

    /**
     * @brief Remove an indentation from print_text output.
     */
    void unindent();
}