#pragma once

// includes
#include <string>
#include <iostream>

using std::string, std::cout, std::endl;

/** Print a red error message to the terminal. 
 * @param str the string to be printed. 
 */
inline void print_err(string str) {
    cout << "\033[1;31m" << str << "\033[0m\n" << endl;
}