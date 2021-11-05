#pragma once

// includes
#include <string>
#include <vector>

// my own includes
#include "../Atom.cpp"

using std::vector, std::string, std::cout, std::endl;

class Writer {
public: 
    std::ofstream file;

    /** Constructor for the Writer class. 
     * @param filename the name of the output file
     */
    Writer(string filename) {
        this->filename = filename;
        file.open(filename, std::ios::trunc);

        // check if file was succesfully opened
        if (!file.is_open()) {
            print_err("Could not open file \"" + filename + "\"");
            exit(1);
        }
    };

    Writer(){};

    virtual void write(vector<Atom*>* protein_atoms, vector<Atom*>* hydration_atoms) {}

    void close() {file.close();}

    string get_filename() {return filename;}

private:
    string filename;
};