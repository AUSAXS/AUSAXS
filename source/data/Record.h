#pragma once

// includes
#include <string>
#include "boost/format.hpp"

#include "Tools.cpp"

using std::string, boost::format;

class Record {
public: 
    enum RecordType {HEADER, ATOM, TERMINATE, FOOTER, NOTYPE};
    
    virtual void parse_pdb(const string s) = 0;
    virtual RecordType get_type() const = 0;
    virtual string as_pdb() const = 0;

    static RecordType get_type(string s) {
        if (type_map.count(s) == 1) {
            return type_map.at(s);
        }
        print_err((format("Error in Record::get_type: Could not determine type \"%1%\"") % s).str());
        exit(1);
    }

private:
    static const inline std::map<string, Record::RecordType> type_map = {
        {"ATOM  ", ATOM}, {"HETATM", ATOM},
        {"TER   ", TERMINATE}, 
        {"HEADER", HEADER}, {"TITLE ", HEADER}, {"COMPND", HEADER}, {"SOURCE", HEADER}, {"KEYWDS", HEADER}, 
        {"EXPDTA", HEADER}, {"AUTHOR", HEADER}, {"REVDAT", HEADER}, {"JRNL  ", HEADER}, {"REMARK", HEADER}, 
        {"DBREF ", HEADER}, {"SEQRES", HEADER}, {"FORMUL", HEADER}, {"HELIX ", HEADER}, {"SHEET ", HEADER}, 
        {"SSBOND", HEADER}, {"CRYST1", HEADER}, {"ORIGX1", HEADER}, {"ORIGX2", HEADER}, {"ORIGX3", HEADER}, 
        {"SCALE1", HEADER}, {"SCALE2", HEADER}, {"SCALE3", HEADER}, 
        {"CONECT", FOOTER}, {"MASTER", FOOTER}, {"END   ", FOOTER}};
};