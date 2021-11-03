#pragma once

// includes
#include <string>
#include <vector>

// my own includes
#include "Atom.cpp"

class Reader {
public: 
    std::ifstream file;

    /** Constructor for the Reader class. 
     * @param filename the name of the input file
     */
    Reader(std::string filename) {
        this->filename = filename;
        loc = 0;
        file.open(filename);

        // check if file was succesfully opened
        if (!file.is_open()) {
            perror(("Could not open file \"" + filename + "\"").c_str());
            exit(EXIT_FAILURE);
        }
    };

    Reader(){};

    virtual std::vector<Atom*> read() {return std::vector<Atom*>();};

private:
    std::string filename;
    int loc; // line number
    virtual Atom read_line() {return Atom();};
};