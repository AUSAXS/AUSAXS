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

    void add(string s) {contents += s;}

    string as_pdb() const {return "";}
    string as_pdbml() const {return "";}

private: 
    string contents;
};