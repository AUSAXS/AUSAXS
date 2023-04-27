#include <utility/Console.h>
#include <utility/ConsoleColor.h>

void console::print_warning(std::string_view text) {
    console::print(text, console::color::red);
}

void console::print_success(std::string_view text) {
    console::print(text, console::color::green);
}

void console::print_info(std::string_view text) {
    console::print(text, console::color::lightblue);
}