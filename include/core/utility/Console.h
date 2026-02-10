// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/ConsoleColor.h>

#include <string_view>

namespace ausaxs::console {
    /**
     * @brief Print a critical message. This will always be printed. 
     */
    void print_critical(std::string_view text, color::color text_color = color::red);

    /**
     * @brief Print further critical text. This will always be printed. 
     */
    void print_text_critical(std::string_view text, color::color text_color = color::white);

    /**
     * @brief Print a warning message. This will only be printed if warnings are enabled. 
     */
    void print_warning(std::string_view text, color::color text_color = color::red);

    /**
     * @brief Print a success message. This will only be printed if verbose output is enabled.
     */
    void print_success(std::string_view text, color::color text_color = color::green);

    /**
     * @brief Print a failure message. This will only be printed if verbose output is enabled.
     */
    void print_failure(std::string_view text, color::color text_color = color::red);

    /**
     * @brief Print a info message. This will only be printed if verbose output is enabled.
     */
    void print_info(std::string_view text, color::color text_color = color::lightblue);

    /**
     * @brief Print a text message. This will only be printed if verbose output is enabled.
     *        Automatic indentation based on the current indentation level will be added to the text.
     */
    void print_text(std::string_view text, color::color text_color = color::white);

    /**
     * @brief Print a minor text message. This will only be printed if minor messages are enabled.
     *        Automatic indentation based on the current indentation level will be added to the text.
     */
    void print_text_minor(std::string_view text, color::color text_color = color::white);

    /**
     * @brief Add another indentation to print_text output.
     */
    void indent(int level = 1);

    /**
     * @brief Remove an indentation from print_text output.
     */
    void unindent(int level = 1);
}