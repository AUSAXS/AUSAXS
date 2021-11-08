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
    vector<shared_ptr<Atom>> get_atoms() {
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
    void add(shared_ptr<Atom> r) {
        contents.push_back(r);
    }

    /**
     * @brief Add a Terminate record to this File. 
     * @param r Terminate to be added. 
     */
    void add(shared_ptr<Terminate> r) {
        contents.push_back(r);
    }

    /**
     * @brief Add a header or footer record to this File. 
     * @param type HEADER or FOOTER.
     * @param s string to be added. 
     */
    void add(string type, string s) {
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
    string as_pdb() {
        string s;
        for (int i = 0; i < contents.size(); i++) {
            s += contents[i]->as_pdb();
        }
        return s;
    }

    /**
     * @brief Create a string .xml representation of this File. 
     * @return The .xml string representation. 
     */
    string as_pdbml() {
        string s;
        for (int i = 0; i < contents.size(); i++) {
            s += contents[i]->as_xml();
        }
        return s;
    }

    static Record::RecordType get_type(string s) {
        if (type_map.count(s) == 1) {
            return type_map.at(s);
        }
        print_err((format("Error in File::get_type: Could not determine type \"%1%\"") % s).str());
        exit(1);
    }

    static const inline std::map<string, Record::RecordType> type_map = {
        {"ATOM  ", Record::ATOM}, {"HETATM", Record::ATOM},
        {"TER   ", Record::TERMINATE}, 
        {"HEADER", Record::HEADER}, {"TITLE", Record::HEADER}, {"COMPND", Record::HEADER}, {"SOURCE", Record::HEADER}, {"KEYWDS", Record::HEADER}, 
        {"EXPDTA", Record::HEADER}, {"AUTHOR", Record::HEADER}, {"REVDAT", Record::HEADER}, {"JRNL  ", Record::HEADER}, {"REMARK", Record::HEADER}, 
        {"DBREF ", Record::HEADER}, {"SEQRES", Record::HEADER}, {"FORMUL", Record::HEADER}, {"HELIX ", Record::HEADER}, {"SHEET ", Record::HEADER}, 
        {"SSBOND", Record::HEADER}, {"CRYST1", Record::HEADER}, {"ORIGX1", Record::HEADER}, {"ORIGX2", Record::HEADER}, {"ORIGX3", Record::HEADER}, 
        {"SCALE1", Record::HEADER}, {"SCALE2", Record::HEADER}, {"SCALE3", Record::HEADER}, 
        {"CONECT", Record::FOOTER}, {"MASTER", Record::FOOTER}, {"END   ", Record::FOOTER}};
};