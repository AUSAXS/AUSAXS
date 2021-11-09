#pragma once

// includes
#include <string>
#include <vector>
#include <utility>

// my own includes
#include "data/Record.h"
#include "data/Terminate.cpp"
#include "data/Header.cpp"
#include "data/Footer.cpp"
#include "data/Atom.cpp"

using std::vector, std::string, std::cout, std::endl, std::unique_ptr, std::shared_ptr; 

class File {
public: 
    vector<shared_ptr<Record>> contents;
    Header header;
    Footer footer;

    /** 
     * @brief Constructor for the File class. 
     */
    File() {};

    /**
     * @brief Get the atoms contained in this File. 
     * @return The atoms from this file. 
     */
    vector<shared_ptr<Atom>> get_atoms() const {
        vector<shared_ptr<Atom>> atoms;
        for (auto const& r : contents) {
            if (r->get_type() == Record::ATOM) {
                atoms.push_back(std::static_pointer_cast<Atom>(r));
            }
        }
        return atoms;
    }

    /** 
     * @brief Add an Atom record to this File. 
     * @param r Atom to be added. 
     */
    void add(const shared_ptr<Atom> r) {
        contents.push_back(r);
    }

    /**
     * @brief Add a Terminate record to this File. 
     * @param r Terminate to be added. 
     */
    void add(const shared_ptr<Terminate> r) {
        contents.push_back(r);
    }

    /**
     * @brief Add a header or footer record to this File. 
     * @param type HEADER or FOOTER.
     * @param s string to be added. 
     */
    void add(const string type, const string s) {
        if (type == "HEADER") {
            header.add(s);
        } else if (type == "FOOTER") {
            footer.add(s);
        } else {
            print_err("Error in File::add: string " + type + " is not \"HEADER\" or \"FOOTER\"!");
            exit(1);
        }
    }

    /**
     * @brief Create a string .pdb representation of this File.
     * @return The .pdb string representation. 
     */
    string as_pdb() const {
        string s;
        s += header.get();
        for (int i = 0; i < contents.size(); i++) {
            s += contents[i]->as_pdb();
        }
        s += footer.get();
        return s;
    }

    /**
     * @brief Create a string .xml representation of this File. 
     * @return The .xml string representation. 
     */
    string as_pdbml() const {
        print_err("Error in File::as_pdbml: Not implemented.");
        exit(1);
    }
};