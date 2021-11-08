#pragma once

// includes
#include <string>
#include <vector>
#include "data/Record.h"

using std::string, std::vector;

class Footer : Record {
public: 
    Footer(){}

    void parse_pdb(string s) {}
    void parse_xml(string s) {}

    RecordType get_type() override {return FOOTER;}

    string as_pdb() const {return "";}
    string as_pdbml() const {return "";}

    void add(string s) {contents += s;}

private: 
    string contents;
};