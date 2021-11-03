#pragma once

// includes
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>

template <typename T>
void save(std::vector<T> v, std::string path) {
    std::filesystem::path p(path);
    std::filesystem::create_directories(p.remove_filename().string()); // create the directories if they do not exist
    std::ofstream file(p.string());
    if (!file.is_open()) {
        perror(("Could not create file at " + p.string()).c_str());
        exit(EXIT_FAILURE);
    }
    for (T e : v) {
        file << e << std::endl;
    }
    file.close();
}