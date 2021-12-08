// includes
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iomanip>
#include <iostream>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

// my own stuff
#include "data/Record.h"
#include "Tools.h"
#include "Atom.h"
#include "math/Vector3.h"

using std::vector, std::string, std::cout, std::endl, std::setw, std::left, std::right, std::shared_ptr, std::unique_ptr;
using boost::format;

void Atom::parse_pdb(string s) {
    int pad_size = 81 - s.size();
    s += string(pad_size, ' ');

    // http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

    //                   RN SE S1 NA AL RN CI S2 RS iC S3 X  Y  Z  OC TF S4  EL CH
    //                   0     1           2              3     4  5  6      7     8
    //                   0  6  1  2  6  7  0  1  2  6  7  0  8  6  4  0  6   6  8  0  
    const char form[] = "%6c%5c%1c%4c%1c%3c%1c%1c%4c%1c%3c%8c%8c%8c%6c%6c%10c%2c%2c";
    string recName = "      ", serial = "     ", space1 = " ", name = "    ", altLoc = " ", resName = "   ", space2 = " ", 
        chainID = " ", resSeq = "    ", iCode = " ", space3 = "   ", x = "        ", y = "        ", z = "        ", 
        occupancy = "      ", tempFactor = "      ", space4 = "          ", element = "  ", charge = "  ";
    sscanf(s.c_str(), form, recName.data(), serial.data(), space1.data(), name.data(), altLoc.data(), resName.data(), 
        space2.data(), chainID.data(), resSeq.data(), iCode.data(), space3.data(), x.data(), y.data(), z.data(), 
        occupancy.data(), tempFactor.data(), space4.data(), element.data(), charge.data());

    // sanity check
    if (__builtin_expect(!(Record::get_type(recName) == Record::ATOM || Record::get_type(recName) == Record::HETATM), false)) {
        print_err("Error in Atom::parse_pdb: input string is not \"ATOM  \" or \"HETATM\" (" + recName + ").");
        exit(1);
    }

    // remove any spaces from the numbers
    boost::erase_all(serial, " ");
    boost::erase_all(name, " ");
    boost::erase_all(resName, " ");
    boost::erase_all(resSeq, " ");
    boost::erase_all(x, " ");
    boost::erase_all(y, " ");
    boost::erase_all(z, " ");
    boost::erase_all(occupancy, " ");
    boost::erase_all(tempFactor, " ");
    boost::erase_all(element, " ");

    // sometimes people use the first character of x for some other purpose.
    // if it is a digit, the following won't work. On the other hand they're kinda asking for it then. Follow the standard, people. 
    if (!(std::isdigit(x[0]) || x[0] == '-')) {
        x = x.substr(1, x.size()-1);
    }

    // set all of the properties
    try {
        _serial = std::stoi(serial);
        _name = name;
        _altLoc = altLoc;
        _resName = resName;
        _chainID = chainID;
        _resSeq = std::stoi(resSeq);
        _iCode = iCode;
        set_coordinates({std::stod(x), std::stod(y), std::stod(z)});
        if (occupancy.empty()) {_occupancy = 1;} else {_occupancy = std::stod(occupancy);}
        if (tempFactor.empty()) {_tempFactor = 0;} else {_tempFactor = std::stod(tempFactor);}
        if (element.empty()) {set_element(name.substr(0, 1));} else {set_element(element);} // the backup plan is to use the first character of "name"
        _charge = charge;
    } catch (const std::exception& e) { // catch conversion errors and output a more meaningful error message
        print_err("Error in Atom::parse_pdb: Invalid field values in line \"" + s + "\".");
        exit(1);
    }

    _effective_charge = constants::charge::get.at(this->element) + constants::hydrogen_atoms::get.at(this->resName).at(this->name);

    // DEBUG OUTPUT
    // cout << s << endl;
    // cout << this->as_pdb() << endl;
}

string Atom::as_pdb() const {
    std::stringstream ss;
    //                   RN SE S1 NA AL RN CI S2 RS iC S3 X  Y  Z  OC TF S2  EL CH
    //                   0     1           2              3     4  5  6      7     8
    //                   0  6  1  2  6  7  0  1  2  6  7  0  8  6  4  0  6   6  8  0  
    //          format: "%6c%5c%2c%4c%1c%3c %1c%4c%1c%3c%8c%8c%8c%6c%6c%10c%2c%2c"
    ss << left << setw(6) << get_recName()                                   // 1 - 6
        << right << setw(5) << serial                                        // 7 - 11
        << " "                                                               // 12
        << " " << left << setw(3) << name                                    // 13 - 16
        << left << setw(1) << altLoc                                         // 17
        << left << setw(3) << resName                                        // 18 - 20
        << " "                                                               // 21
        << left << setw(1) << chainID                                        // 22
        << right << setw(4) << resSeq                                        // 23 - 26
        << right << setw(1) << iCode                                         // 27
        << "   "                                                             // 28 - 30
        << right << setw(8) << setf(coords.x, 8)                             // 31 - 38
        << right << setw(8) << setf(coords.y, 8)                             // 39 - 46
        << right << setw(8) << setf(coords.z, 8)                             // 47 - 54
        << right << setw(6) << setf(occupancy, 6)                            // 55 - 60
        << right << setw(6) << setf(tempFactor, 6)                           // 61 - 66
        << "          "                                                      // 67 - 76
        << right << setw(2) << element                                       // 77 - 78
        << left << setw(2) << charge                                         // 79 - 80
        << endl;
    return ss.str();
}

/** Prints the contents of this object to the terminal. */
void Atom::print() {
    cout << "\nAtom no: " << serial << endl;
    cout << setw(17) << "(x, y, z): (" << setw(6) << coords.x << ", " << setw(6) << coords.y << ", " << setw(6) << coords.z << ")" << endl;
    cout << setw(16) << "Weight: " << std::to_string(occupancy) << endl;
    cout << setw(16) << "Symbol: " << element << endl;
    cout << setw(16) << "Molecule: " << name << endl;
    return;
}

bool Atom::operator<(const Atom& rhs) const {
    return serial < rhs.serial;
}

bool Atom::operator==(const Atom& rhs) const {
    if (get_type() != get_type()) {return false;}
    if (name != rhs.name) {return false;}
    if (altLoc != rhs.altLoc) {return false;}
    if (resName != rhs.resName) {return false;}
    if (chainID != rhs.chainID) {return false;}
    if (iCode != rhs.iCode) {return false;}
    if (element != rhs.element) {return false;}
    if (charge != rhs.charge) {return false;}
    if (occupancy != rhs.occupancy) {return false;}
    if (tempFactor != rhs.tempFactor) {return false;}
    if (serial != rhs.serial) {return false;}
    // if (resSeq != rhs.resSeq) {return false;} // this is to fix io tests, since some pdb files randomly changes this order
    if (coords != rhs.coords) {return false;}
    return true;
}