#pragma once

// includes
#include <string>
#include <iostream>

/** Print a red error message to the terminal. 
 * @param str the string to be printed. 
 */
void print_err(std::string str) {
    std::cout << "\033[1;31mbold" << str << "\033[0m\n" << std::endl;
}