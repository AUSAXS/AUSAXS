#pragma once

// includes
#include <string>
#include <vector>

// my own includes
#include "data/Atom.cpp"

using std::vector, std::string, std::cout, std::endl;

class Writer {
public: 
    std::ofstream output;

    /** Constructor for the Writer class. 
     * @param filename the name of the output file
     */
    Writer(string filename) {
        this->filename = filename;
        output.open(filename, std::ios::trunc);

        // check if file was succesfully opened
        if (!output.is_open()) {
            print_err("Could not open file \"" + filename + "\"");
            exit(1);
        }
    };

    Writer(){};

    virtual void write(shared_ptr<File> file) {}

    void close() {output.close();}

    string get_filename() {return filename;}

protected:
    string filename;
};