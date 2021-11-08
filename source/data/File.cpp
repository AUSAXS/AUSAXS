#pragma once

// includes
#include <string>
#include <vector>

// my own includes
#include "data/Record.h"
#include "data/Terminate.cpp"
#include "data/Header.cpp"
#include "data/Footer.cpp"
#include "Atom.cpp"

using std::vector, std::string, std::cout, std::endl;

class File {
public: 
    /** 
     * @brief Constructor for the File class. 
     */
    File() {};

    /** 
     * @brief Add an Atom record to this File. 
     * @param r Atom to be added. 
     */
    void add(Atom* r) {
        contents.push_back(r);
    }

    /**
     * @brief Add a Terminate record to this File. 
     * @param r Terminate to be added. 
     */
    void add(Terminate* r) {
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
            print_err("ERROR: string " + type + " is not \"HEADER\" or \"FOOTER\"!");
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

    vector<Record*> contents;
    Header header;
    Footer footer;

    enum Type {HEADER, ATOM, TERMINATE, FOOTER};
    static Type get_type(string s) {
        if (type_map.count(s) == 1) {
            return type_map.at(s);
        }
        print_err((format("ERROR: Could not determine type \"%1%\"") % s).str());
        exit(1);
    }

    static const inline std::map<string, Type> type_map = {{"ATOM  ", ATOM}, {"HETATM", ATOM}, {"TER   ", TERMINATE}, 
        {"HEADER", HEADER}, {"TITLE", HEADER}, {"COMPND", HEADER}, {"SOURCE", HEADER}, {"KEYWDS", HEADER}, {"EXPDTA", HEADER},
        {"AUTHOR", HEADER}, {"REVDAT", HEADER}, {"JRNL  ", HEADER}, {"REMARK", HEADER}, {"DBREF ", HEADER}, {"SEQRES", HEADER},
        {"FORMUL", HEADER}, {"HELIX ", HEADER}, {"SHEET ", HEADER}, {"SSBOND", HEADER}, {"CRYST1", HEADER}, {"ORIGX1", HEADER},
        {"ORIGX2", HEADER}, {"ORIGX3", HEADER}, {"SCALE1", HEADER}, {"SCALE2", HEADER}, {"SCALE3", HEADER}, 
        {"CONECT", FOOTER}, {"MASTER", FOOTER}, {"END   ", FOOTER}};
};