#pragma once

// includes
#include <string>
#include <vector>

// my own includes
#include "data/File.cpp"

using std::vector, std::string, std::cout, std::endl, std::unique_ptr, std::shared_ptr;

class Reader {
public: 
    /** Constructor for the Reader class. 
     * @param filename the name of the input file
     */
    Reader(string filename) {
        this->filename = filename;
        input.open(filename);

        // check if file was succesfully opened
        if (!input.is_open()) {
            print_err("Could not open file \"" + filename + "\"");
            exit(1);
        }
    };

    Reader(){}

    virtual unique_ptr<File> read() {return std::make_unique<File>();}

    void close() {input.close();}

    string get_filename() {return filename;}

protected:
    std::ifstream input;
    string filename;
};