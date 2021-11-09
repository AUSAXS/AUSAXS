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
    File(string filename) {this->filename = filename;}

    Header header;
    Footer footer;
    vector<shared_ptr<Record>> contents;

    /**
     * @brief Update the contents of this File to reflect the input Protein.
     * @param protein the protein to use.
     */
    virtual void update(vector<shared_ptr<Atom>> protein_atoms, vector<shared_ptr<Atom>> hydration_atoms) = 0;

    /**
     * @brief Read the file backing this File object. 
     */
    virtual void read() = 0;

    /**
     * @brief write this File to disk. 
     * @param path the output path.
     */
    virtual void write(const string path) const = 0;

    /**
     * @brief Get the protein atoms contained in this File. 
     * @return The atoms from this file. 
     */
    virtual vector<shared_ptr<Atom>> get_protein_atoms() const {
        vector<shared_ptr<Atom>> atoms;
        for (auto const& r : contents) {
            if (r->get_type() == Record::ATOM) {
                shared_ptr<Atom> a = std::static_pointer_cast<Atom>(r);
                if (!a->is_water()) {
                    atoms.push_back(a);
                }
            }
        }
        return atoms;
    };

    /**
     * @brief Get the hydration atoms contained in this File. 
     * @return The atoms from this file. 
     */
    virtual vector<shared_ptr<Atom>> get_hydration_atoms() const {
        vector<shared_ptr<Atom>> atoms;
        for (auto const& r : contents) {
            if (r->get_type() == Record::ATOM) {
                shared_ptr<Atom> a = std::static_pointer_cast<Atom>(r);
                if (a->is_water()) {
                    atoms.push_back(a);
                }
            }
        }
        return atoms;
    };

    /**
     * @brief Get the atoms contained in this File. 
     * @return A pair of [protein, hydration] atoms contained in this File. 
     */
    std::pair<vector<shared_ptr<Atom>>, vector<shared_ptr<Atom>>> get_atoms() const {
        vector<shared_ptr<Atom>> protein_atoms;
        vector<shared_ptr<Atom>> hydration_atoms;
        for (auto const& r : contents) {
            if (r->get_type() == Record::ATOM) {
                shared_ptr<Atom> a = std::static_pointer_cast<Atom>(r);
                if (a->is_water()) {
                    hydration_atoms.push_back(a);
                } else {
                    protein_atoms.push_back(a);
                }
            }
        }
        return std::make_pair(protein_atoms, hydration_atoms);
    }

    /**
     * @brief Create a string representation of this File.
     * @return The string representation. 
     */
    string as_pdb() const {
        string s;
        s += header.get();
        for (int i = 0; i < contents.size(); i++) {
            s += contents[i]->as_pdb();
        }
        s += footer.get();
        return s;
    };

    /** 
     * @brief Add an Atom record to this File. 
     * @param r Atom to be added. 
     */
    virtual void add(const shared_ptr<Atom> r) {
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

protected:
    string filename;
};