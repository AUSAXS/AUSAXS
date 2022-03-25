#pragma once

#include <string>
#include <sstream>
#include <map>
#include <math.h>

#include <Exceptions.h>

class Record {
  public: 
    enum RecordType {HEADER, ATOM, HETATM, TERMINATE, FOOTER, NOTYPE};
    virtual ~Record() {}
    
    virtual void parse_pdb(const string s) = 0;
    virtual RecordType get_type() const = 0;
    virtual string as_pdb() const = 0;

    static RecordType get_type(string s) {
        if (type_map.count(s) == 1) {
            return type_map.at(s);
        }
        throw except::parse_error("Error in Record::get_type: Could not determine type \"" + s + "\"");
    }

  private:
    inline static const std::map<string, RecordType> type_map = {
        {"ATOM"  , ATOM}, {"HETATM", HETATM},
        {"TER"   , TERMINATE}, 
        {"HEADER", HEADER}, {"TITLE" , HEADER}, {"COMPND", HEADER}, {"SOURCE", HEADER}, {"KEYWDS", HEADER}, 
        {"EXPDTA", HEADER}, {"AUTHOR", HEADER}, {"REVDAT", HEADER}, {"JRNL"  , HEADER}, {"REMARK", HEADER}, 
        {"DBREF" , HEADER}, {"SEQRES", HEADER}, {"FORMUL", HEADER}, {"HELIX" , HEADER}, {"SHEET" , HEADER}, 
        {"SSBOND", HEADER}, {"CRYST1", HEADER}, {"ORIGX1", HEADER}, {"ORIGX2", HEADER}, {"ORIGX3", HEADER}, 
        {"SCALE1", HEADER}, {"SCALE2", HEADER}, {"SCALE3", HEADER}, {"HET"   , HEADER}, {"HETNAM", HEADER},
        {"HETSYN", HEADER}, {"FORMUL", HEADER}, {"CISPEP", HEADER}, {"SITE"  , HEADER},
        {"CONECT", FOOTER}, {"MASTER", FOOTER}, {"END"   , FOOTER}};
};

// Fixed-length printing of numbers. std::setprecision does *not* count leading zeros, which breaks our strict formatting.
struct __setp {
    double number;
    int prec;
};

inline std::ostream& operator<<(std::ostream& os, const __setp& obj) {
    os.precision(obj.prec);
    os << obj.number;
    return os;
}

/**
 * @brief Print a fixed-width number in a given stream.
 * @param number number to be printed. 
 * @param p total width including signs and decimal separator
 */
inline __setp setf(const double number, const int precision) {
    // IF THIS EVER BREAKS AGAIN, CONVERT IT TO A SIMPLE LOOP OVER string(number) INSTEAD. COULD EVEN INCLUDE MINWIDTH AS WELL
    __setp setter;
    setter.number = number;
    int p = 1; // the dot will always take up a slot

    string num = std::to_string(number);
    if (num[0] == '-') {
        p++; // -1 since the sign takes up a slot
    }

    if (std::floor(std::abs(number)) != 0) { // if the number is not ±0.xxx
        setter.prec = precision - p;
        return setter;
    }

    p++; // another slot is taken by 0
    size_t fsig = num.find_first_not_of("0", p); // find first significant decimal after the dot
    if (fsig != string::npos) { // number is of the form 0.00312...
        p = fsig;
    }
    setter.prec = precision - p;
    return setter;
}
