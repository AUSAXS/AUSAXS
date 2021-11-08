#pragma once

// includes
#include <string>

using std::string;

class Record {
public: 
    virtual void parse_pdb(const string s) {}
    virtual void parse_xml(const string s) {}

    string as_pdb() const {return "";}
    string as_xml() const {return "";}
};