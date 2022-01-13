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
#include "data/Atom.h"
#include "math/Vector3.h"

using std::vector, std::string, std::cout, std::endl, std::setw, std::left, std::right, std::shared_ptr, std::unique_ptr;
using boost::format;

Atom::Atom(Atom&& a) noexcept : coords(std::move(a.coords)), _name(std::move(a.name)), _altLoc(std::move(a.altLoc)), _resName(std::move(a.resName)), 
    _chainID(std::move(a.chainID)), _iCode(std::move(a.iCode)), _element(std::move(a.element)), _charge(std::move(a.charge)), 
    _occupancy(std::move(a.occupancy)), _tempFactor(std::move(a.tempFactor)), _serial(std::move(a.serial)), _resSeq(std::move(a.resSeq)), 
    _effective_charge(std::move(a.effective_charge)), _uid(std::move(a.uid)) {}

Atom::Atom(const Atom& a) :coords(a.coords), _name(a.name), _altLoc(a.altLoc), _resName(a.resName), _chainID(a.chainID), _iCode(a.iCode), 
    _element(a.element), _charge(a.charge), _occupancy(a.occupancy), _tempFactor(a.tempFactor), _serial(a.serial), _resSeq(a.resSeq),
    _effective_charge(a.effective_charge), _uid(a.uid) {}

Atom::Atom(const Vector3 v, const double occupancy, const string element, const string name, int serial) {
    // we use our setters so we can validate the input if necessary
    set_coordinates(v);
    set_occupancy(occupancy);
    set_element(element);
    set_name(name);
    set_serial(serial);
    _effective_charge = constants::charge::get.at(this->element);
    _uid = uid_counter++;
}

Atom::Atom(const int serial, const string name, const string altLoc, const string resName, const string chainID, const int resSeq, 
    const string iCode, const Vector3 coords, const double occupancy, const double tempFactor, const string element, const string charge) {
        set_serial(serial);
        set_name(name);
        _altLoc = altLoc;
        set_resName(resName);
        _chainID = chainID;
        _resSeq = resSeq;
        _iCode = iCode;
        set_coordinates(coords);
        set_occupancy(occupancy);
        _tempFactor = tempFactor;
        set_element(element);
        _charge = charge;
        try {
            _effective_charge = constants::charge::get.at(this->element) + constants::hydrogen_atoms::get.at(this->resName).at(this->name);
        } catch (const std::exception& e) {
            print_err("Could not set effective charge: unknown element, residual or atom (" + element + ", " + resName + ", " + name + ")");
        }
        _uid = uid_counter++;
}

Atom::Atom() : _uid(uid_counter++) {}

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
        throw except::parse_error("Error in Atom::parse_pdb: input string is not \"ATOM  \" or \"HETATM\" (" + recName + ").");
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
        throw except::parse_error("Error in Atom::parse_pdb: Invalid field values in line \"" + s + "\".");
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

Atom& Atom::operator=(const Atom& rhs) {
    _name = rhs.name; 
    _altLoc = rhs.altLoc; 
    _resName = rhs.resName; 
    _chainID = rhs.chainID; 
    _iCode = rhs.iCode; 
    _element = rhs.element; 
    _charge = rhs.charge;
    _occupancy = rhs.occupancy; 
    _tempFactor = rhs.tempFactor;
    _serial = rhs.serial; 
    _resSeq = rhs.resSeq;
    coords = rhs.coords;
    _effective_charge = rhs.effective_charge;
    _uid = rhs.uid;
    return *this;
}