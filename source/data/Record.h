#pragma once

// includes
#include <string>

using std::string;


class Record {
public: 
    enum RecordType {HEADER, ATOM, TERMINATE, FOOTER, NOTYPE};
    
    virtual void parse_pdb(const string s) {}
    virtual void parse_xml(const string s) {}
    virtual RecordType get_type() {return RecordType::NOTYPE;}
    virtual string as_pdb() const {return "";}
    virtual string as_xml() const {return "";}
};