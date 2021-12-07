#pragma once

// includes
#include <string>
#include <vector>
#include <utility>
#include <fstream>

// my own includes
#include "data/Record.h"
#include "data/Terminate.h"
#include "data/Header.h"
#include "data/Footer.h"
#include "data/Atom.h"
#include "data/Hetatom.h"
#include "data/File.h"
#include "settings.h"

using std::vector, std::string, std::cout, std::endl, std::unique_ptr, std::shared_ptr; 

class PDB_file : public File {
public: 
    /** 
     * @brief Constructor for the PDB_file class. 
     */
    PDB_file(string filename) : File(filename) {
        read();
    }
    ~PDB_file() override {}
    
    /**
     * @brief write this File to disk. 
     * @param path the output path.
     */
    void write(string path) const override {
        std::ofstream output(path);
        if (!output.is_open()) {throw std::ios_base::failure("Error in PDB_file::write: Could not open file \"" + path + "\"");}
        output << as_pdb() << std::flush;
        output.close();
        cout << "Output written to file " + path + "." << endl;
    }

private:
    /**
     * @brief Read the file backing this File object. 
     */
    void read() override {
        // check if file was succesfully opened
        std::ifstream input(filename);
        if (!input.is_open()) {throw std::ios_base::failure("Error in PDB_file::read: Could not open file \"" + filename + "\"");}

        string line; // placeholder for the current line
        while(getline(input, line)) {
            string type = line.substr(0, 6); // read the first 6 characters
            switch(Record::get_type(type)) {
                case Record::RecordType::HETATM: {
                    Hetatom atom;
                    atom.parse_pdb(line);
                    add(atom);
                    break;
                } case Record::RecordType::ATOM: {
                    Atom atom;
                    atom.parse_pdb(line);
                    add(atom);
                    break;
                } case Record::RecordType::TERMINATE: {
                    Terminate term;
                    term.parse_pdb(line);
                    add(term);
                    break;
                } case Record::RecordType::HEADER: {
                    add("HEADER", line);
                    break;
                } case Record::RecordType::FOOTER: {
                    add("FOOTER", line);
                    break;
                } default: {
                    throw std::ios_base::failure("Error in PDB_file::read: Malformed input file - unrecognized type.");
                }
            };
        }
        input.close();
    }
};