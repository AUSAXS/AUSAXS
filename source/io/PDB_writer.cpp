#pragma once

// includes
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// my own includes
#include "Writer.h"
#include "data/File.cpp"

using std::setw, std::left, std::right;

class PDB_writer : public Writer {
public: 
    /** Constructor for the PDBML_writer class. 
     * @param filename the name of the input file
     */
    PDB_writer(string filename) : Writer(filename) {};

    void write(shared_ptr<File> file) override {output << file->as_pdb() << std::flush;}
};