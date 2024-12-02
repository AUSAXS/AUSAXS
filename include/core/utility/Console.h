#pragma once

#include <string_view>

namespace ausaxs::console {
    /**
     * @brief Enable logging. This will replicate all console output to a log file.
     */
    void enable_logging();

    /**
     * @brief Disable logging.
     */
    void disable_logging();

    /**
     * @brief Print a critical message. The text will be red in the terminal. 
     *        This will always be printed. 
     */
    void print_critical(std::string_view text);

    /**
     * @brief Print further critical text. The text will be white in the terminal.
     *        This will always be printed. 
     */
    void print_text_critical(std::string_view text);

    /**
     * @brief Print a warning message. The text will be red in the terminal. 
     *        This will only be printed if verbose output is enabled. 
     */
    void print_warning(std::string_view text);

    /**
     * @brief Print a success message. The text will be green in the terminal. 
     *        This will only be printed if verbose output is enabled.
     */
    void print_success(std::string_view text);

    /**
     * @brief Print a info message. The text will be blue in the terminal. 
     *        This will only be printed if verbose output is enabled.
     */
    void print_info(std::string_view text);

    /**
     * @brief Print a text message. The text will be white in the terminal. 
     *        This will only be printed if verbose output is enabled.
     *        Automatic indentation based on the current indentation level will be added to the text.
     */
    void print_text(std::string_view text);

    /**
     * @brief Print a minor text message. The text will be white in the terminal. 
     *        This will only be printed if minor messages are enabled.
     *        Automatic indentation based on the current indentation level will be added to the text.
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