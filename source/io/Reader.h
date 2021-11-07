#pragma once

// includes
#include <string>
#include <vector>

// my own includes
#include "Atom.cpp"

using std::vector, std::string, std::cout, std::endl;

class Reader {
public: 
    std::ifstream file;

    /** Constructor for the Reader class. 
     * @param filename the name of the input file
     */
    Reader(string filename) {
        this->filename = filename;
        file.open(filename);

        // check if file was succesfully opened
        if (!file.is_open()) {
            print_err("Could not open file \"" + filename + "\"");
            exit(1);
        }
    };

    Reader(){}

    virtual vector<Atom*> read() {return vector<Atom*>();}

    void close() {file.close();}

    string get_filename() {return filename;}

private:
    string filename;
};