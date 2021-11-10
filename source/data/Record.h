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

// Fixed-length printing of numbers. std::setprecision does *not* count leading zeros, which breaks our strict formatting.
// Based on https://stackoverflow.com/a/66691675 
struct __setp {
    double number;
    int prec;
};

std::ostream& operator<<(std::ostream& os, const __setp& obj) {
    os.precision(obj.prec);
    os << obj.number;
    return os;
}

__setp setp(double number, int p) {
    __setp setter;
    setter.number = number;

    string num = std::to_string(number);
    int fdec = num.find("0.")+1; // find first decimal after a 0
    if (fdec == string::npos) { // if it does not start with 0, there are no issues
        setter.prec = p;
        return setter;
    }

    int fsig = num.find_first_not_of("0", fdec); // find first significant decimal
    if (fsig == string::npos) { // number is of the form 0.0000...
        setter.prec = p - fdec; 
    } else { // number is of the form 0.00314...
        setter.prec = p - fsig;
    }
    return setter;
}
