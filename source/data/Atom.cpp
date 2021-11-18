// includes
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iomanip>
#include <iostream>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

// ROOT
#include <TVector3.h>

// my own stuff
#include "data/Record.h"
#include "Tools.h"
#include "Atom.h"

using namespace ROOT;
using std::vector, std::string, std::cout, std::endl, std::setw, std::left, std::right, std::shared_ptr, std::unique_ptr;
using boost::format;

void Atom::parse_pdb(const string s) {
    // http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

    //                   RN SE S1 NA AL RN CI S2 RS iC S3 X  Y  Z  OC TF S2  EL CH
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
    if (!Record::get_type(recName) == Record::ATOM) {
        print_err("Error in Atom::parse_pdb: input string is not \"ATOM  \" or \"HETATM\" (" + recName + ").");
        exit(1);
    }

    // remove any spaces
    boost::erase_all(serial, " ");
    boost::erase_all(resSeq, " ");
    boost::erase_all(x, " ");
    boost::erase_all(y, " ");
    boost::erase_all(z, " ");
    boost::erase_all(occupancy, " ");
    boost::erase_all(tempFactor, " ");

    // sometimes people use the first character of x for some other purpose.
    // if it is a digit, the following won't work. On the other hand they're kinda asking for it then. Follow the standard, people. 
    if (!(std::isdigit(x[0]) || x[0] == '-')) {
        x = x.substr(1, x.size()-1);
    }

    // set all of the properties
    try {
        this->recName = recName;
        this->set_serial(std::stoi(serial));
        this->set_name(name);
        this->altLoc = altLoc;
        this->set_resName(resName);
        this->chainID = chainID;
        this->resSeq = std::stoi(resSeq);
        this->iCode = iCode;
        this->set_coordinates({std::stod(x), std::stod(y), std::stod(z)});
        this->set_occupancy(std::stod(occupancy));
        this->tempFactor = std::stod(tempFactor);
        this->set_element(element);
        this->charge = charge;
    } catch (const std::exception& e) { // catch conversion errors and output a more meaningful error message
        print_err("Error in Atom::parse_pdb: Invalid field values in line \"" + s + "\".");
        exit(1);
    }

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
    ss << left << setw(6) << recName                                         // 1 - 6
        << right << setw(5) << get_serial()                                  // 7 - 11
        << " "                                                               // 12
        << " " << left << setw(3) << name                                    // 13 - 16
        << left << setw(1) << altLoc                                         // 17
        << left << setw(3) << resName                                        // 18 - 20
        << " "                                                               // 21
        << left << setw(1) << chainID                                        // 22
        << right << setw(4) << resSeq                                        // 23 - 26
        << right << setw(1) << iCode                                         // 27
        << "   "                                                             // 28 - 30
        << right << setw(8) << setf(get_x(), 8)                              // 31 - 38
        << right << setw(8) << setf(get_y(), 8)                              // 39 - 46
        << right << setw(8) << setf(get_z(), 8)                              // 47 - 54
        << right << setw(6) << setf(get_occupancy(), 6)                      // 55 - 60
        << right << setw(6) << setf(tempFactor, 6)                           // 61 - 66
        << "          "                                                      // 67 - 76
        << right << setw(2) << get_element()                                 // 77 - 78
        << left << setw(2) << charge                                         // 79 - 80
        << endl;
    return ss.str();
}

/** Prints the contents of this object to the terminal. */
void Atom::print() {
    cout << "\nAtom no: " << serial << endl;
    cout << setw(17) << "(x, y, z): (" << setw(6) << get_x() << ", " << setw(6) << get_y() << ", " << setw(6) << get_z() << ")" << endl;
    cout << setw(16) << "Weight: " << std::to_string(occupancy) << endl;
    cout << setw(16) << "Symbol: " << element << endl;
    cout << setw(16) << "Molecule: " << name << endl;
    return;
}

bool Atom::operator<(const Atom& rhs) const {
    return serial < rhs.get_serial();
}

bool Atom::operator==(const Atom& rhs) const {
    if (recName != rhs.recName) {return false;}
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
    if (resSeq != rhs.resSeq) {return false;}
    if (coords != rhs.coords) {return false;}
    return true;
}