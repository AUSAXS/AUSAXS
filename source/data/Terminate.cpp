#pragma once

// includes
#include <string>
#include "data/Record.h"

using std::string;

class Terminate : public Record {
public: 
    void parse_pdb(string s) {}
    void parse_xml(string s) {}

    RecordType get_type() override {return TERMINATE;}

    string as_pdb() const {return "";}
    string as_pdbml() const {return "";}
};