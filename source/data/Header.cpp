#pragma once

// includes
#include <string>
#include <vector>
#include "data/Record.h"

using std::string, std::vector;

class Header : Record {
public: 
    Header(){}

    void parse_pdb(string s) {}
    void parse_xml(string s) {}

    RecordType get_type() override {return HEADER;}

    string as_pdb() const {return "";}
    string as_pdbml() const {return "";}

    void add(string s) {contents += s;}

private: 
    string contents;
};