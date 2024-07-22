/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <utility/Console.h>
#include <utility/ConsoleColor.h>
#include <settings/GeneralSettings.h>

std::string indentation = "";
void console::indent() {
    indentation += "\t";
}

void console::unindent() {
    if (indentation.empty()) {
        throw std::runtime_error("Cannot unindent console output below 0.");
    }
    indentation.pop_back();
}

void console::print_critical(std::string_view text) {
    console::print(text, console::color::red);
}

void console::print_text_critical(std::string_view text) {
    console::print(std::string(text), console::color::white);
}

void console::print_warning(std::string_view text) {
    if (!settings::general::verbose) {return;}
    console::print(text, console::color::red);
}

void console::print_success(std::string_view text) {
    if (!settings::general::verbose) {return;}
    console::print(text, console::color::green);
}

void console::print_info(std::string_view text) {
    if (!settings::general::verbose) {return;}
    console::print(text, console::color::lightblue);
}

void console::print_text(std::string_view text) {
    if (!settings::general::verbose) {return;}
    console::print(indentation + std::string(text), console::color::white);
}

bool minor_messages = true;
void console::print_text_minor(std::string_view text) {
    if (!minor_messages || !settings::general::verbose) {return;}
    console::print(indentation + std::string(text), console::color::white);
}