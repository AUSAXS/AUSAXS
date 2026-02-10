// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/Console.h>
#include <utility/ConsoleColor.h>
#include <utility/Logging.h>
#include <settings/GeneralSettings.h>

using namespace ausaxs;

std::string indentation = "";
void console::indent(int level) {
    indentation += std::string(level, '\t');
}

void console::unindent(int level) {
    if (indentation.empty()) {
        throw std::runtime_error("Cannot unindent console output below 0.");
    }
    indentation.resize(indentation.size() - level);
}

void console::print_critical(std::string_view text, color::color text_color) {
    logging::log_critical(text);
    console::print(text, text_color);
}

void console::print_text_critical(std::string_view text, color::color text_color) {
    logging::log_critical(text);
    console::print(std::string(text), text_color);
}

void console::print_warning(std::string_view text, color::color text_color) {
    logging::log_console(text);
    if (!settings::general::warnings) {return;}
    console::print(text, text_color);
}

void console::print_success(std::string_view text, color::color text_color) {
    logging::log_console(text);
    if (!settings::general::verbose) {return;}
    console::print(indentation + std::string(text), text_color);
}

void console::print_failure(std::string_view text, color::color text_color) {
    logging::log_console(text);
    if (!settings::general::verbose) {return;}
    console::print(indentation + std::string(text), text_color);
}

void console::print_info(std::string_view text, color::color text_color) {
    logging::log_console(text);
    if (!settings::general::verbose) {return;}
    console::print(text, text_color);
}

void console::print_text(std::string_view text, color::color text_color) {
    logging::log_console(text);
    if (!settings::general::verbose) {return;}
    console::print(indentation + std::string(text), text_color);
}

bool minor_messages = true;
void console::print_text_minor(std::string_view text, color::color text_color) {
    logging::log_console(text);
    if (!minor_messages || !settings::general::verbose) {return;}
    console::print(indentation + std::string(text), text_color);
}