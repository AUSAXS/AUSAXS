#pragma once

// includes
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <boost/format.hpp>

using std::vector, std::string, std::cout, std::endl;

template <typename T>
void save(vector<T> v, string path) {
    std::filesystem::path p(path);
    std::filesystem::create_directories(p.remove_filename().string()); // create the directories if they do not exist
    std::ofstream file(path, std::ios::trunc);
    if (!file.is_open()) {
        perror(("Could not create file at " + path).c_str());
        exit(EXIT_FAILURE);
    }
    for (T e : v) {
        file << e << endl;
    }
    file.close();
}

/** Print a red error message to the terminal. 
 * @param str the string to be printed. 
 */
void print_err(string str) {
    cout << "\033[1;31m" << str << "\033[0m\n" << endl;
}