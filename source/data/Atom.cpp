#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iomanip>
#include <iostream>
#include <boost/algorithm/string.hpp>

#include <data/Record.h>
#include <data/Atom.h>
#include <math/Vector3.h>
#include <constants.h>

using std::vector, std::string, std::cout, std::endl, std::setw, std::left, std::right, std::shared_ptr, std::unique_ptr;

Atom::Atom(const Atom&& a) noexcept : coords(a.coords), name(a.name), altLoc(a.altLoc), resName(a.resName), chainID(a.chainID), iCode(a.iCode), element(a.element), 
    charge(a.charge), occupancy(a.occupancy), tempFactor(a.tempFactor), serial(a.serial), resSeq(a.resSeq), effective_charge(a.effective_charge), uid(a.uid) {}

Atom::Atom(const Atom& a) :coords(a.coords), name(a.name), altLoc(a.altLoc), resName(a.resName), chainID(a.chainID), iCode(a.iCode), element(a.element), 
    charge(a.charge), occupancy(a.occupancy), tempFactor(a.tempFactor), serial(a.serial), resSeq(a.resSeq), effective_charge(a.effective_charge), uid(a.uid) {}

Atom::Atom(const Vector3 v, const double occupancy, const string element, const string name, int serial) {
    // we use our setters so we can validate the input if necessary
    set_coordinates(v);
    set_occupancy(occupancy);
    set_element(element);
    set_name(name);
    set_serial(serial);
    set_effective_charge(constants::charge::get.at(this->element));
    uid = uid_counter++;
}

Atom::Atom(const int serial, const string name, const string altLoc, const string resName, const string chainID, const int resSeq, 
    const string iCode, const Vector3 coords, const double occupancy, const double tempFactor, const string element, const string charge) {
        set_serial(serial);
        set_name(name);
        set_altLoc(altLoc);
        set_resName(resName);
        set_chainID(chainID);
        set_resSeq(resSeq);
        set_iCode(iCode);
        set_coordinates(coords);
        set_occupancy(occupancy);
        set_tempFactor(tempFactor);
        set_element(element);
        set_charge(charge);
        try {
            effective_charge = constants::charge::get.at(this->element) + constants::hydrogen_atoms::get.at(this->resName).at(this->name);
        } catch (const std::exception& e) {
            throw except::invalid_argument("Error in Atom::Atom: Could not set effective charge. Unknown element, residual or atom: (" + element + ", " + resName + ", " + name + ")");
        }
        uid = uid_counter++;
}

Atom::Atom() : uid(uid_counter++) {}

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
        this->serial = std::stoi(serial);
        this->name = name;
        this->altLoc = altLoc;
        this->resName = resName;
        this->chainID = chainID;
        this->resSeq = std::stoi(resSeq);
        this->iCode = iCode;
        set_coordinates({std::stod(x), std::stod(y), std::stod(z)});
        if (occupancy.empty()) {this->occupancy = 1;} else {this->occupancy = std::stod(occupancy);}
        if (tempFactor.empty()) {this->tempFactor = 0;} else {this->tempFactor = std::stod(tempFactor);}
        if (element.empty()) {set_element(name.substr(0, 1));} else {set_element(element);} // the backup plan is to use the first character of "name"
        this->charge = charge;
    } catch (const std::exception& e) { // catch conversion errors and output a more meaningful error message
        throw except::parse_error("Error in Atom::parse_pdb: Invalid field values in line \"" + s + "\".");
    }

    effective_charge = constants::charge::get.at(this->element) + constants::hydrogen_atoms::get.at(this->resName).at(this->name);

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
        << right << setw(8) << setf(coords.x(), 8)                             // 31 - 38
        << right << setw(8) << setf(coords.y(), 8)                             // 39 - 46
        << right << setw(8) << setf(coords.z(), 8)                             // 47 - 54
        << right << setw(6) << setf(occupancy, 6)                            // 55 - 60
        << right << setw(6) << setf(tempFactor, 6)                           // 61 - 66
        << "          "                                                      // 67 - 76
        << right << setw(2) << element                                       // 77 - 78
        << left << setw(2) << charge                                         // 79 - 80
        << endl;
    return ss.str();
}

Record::RecordType Atom::get_type() const {return ATOM;}

double Atom::distance(const Atom& a) const {return coords.distance(a.coords);}
void Atom::translate(const Vector3 v) {coords += v;}
bool Atom::is_water() const {return false;}

void Atom::set_coordinates(const Vector3 v) {coords = v;}
void Atom::set_x(double x) {coords.x() = x;}
void Atom::set_y(double y) {coords.y() = y;}
void Atom::set_z(double z) {coords.z() = z;}
void Atom::set_occupancy(double occupancy) {this->occupancy = occupancy;}
void Atom::set_tempFactor(double tempFactor) {this->tempFactor = tempFactor;}
void Atom::set_altLoc(string altLoc) {this->altLoc = altLoc;}
void Atom::set_serial(int serial) {this->serial = serial;}
void Atom::set_resSeq(int resSeq) {this->resSeq = resSeq;}
void Atom::set_effective_charge(double charge) {effective_charge = charge;}
void Atom::set_chainID(string chainID) {this->chainID = chainID;}
void Atom::set_iCode(string iCode) {this->iCode = iCode;}
void Atom::set_charge(string charge) {this->charge = charge;}
void Atom::set_resName(string resName) {this->resName = resName;}
void Atom::set_name(string name) {this->name = name;}

void Atom::set_element(string element) {
    if (__builtin_expect(constants::mass::atomic.count(element) == 0, false)) { // check that the weight is defined
        throw except::invalid_argument("Error in Atom::set_element: The weight of element " + element + " is not defined.");
    }
    this->element = element;
}

Vector3& Atom::get_coordinates() {return coords;}
const Vector3& Atom::get_coordinates() const {return coords;}
int Atom::get_serial() const {return serial;}
int Atom::get_resSeq() const {return resSeq;}
double Atom::get_occupancy() const {return occupancy;}
double Atom::get_tempFactor() const {return tempFactor;}
double Atom::get_effective_charge() const {return effective_charge;}
string Atom::get_altLoc() const {return altLoc;}
string Atom::get_chainID() const {return chainID;}
string Atom::get_iCode() const {return iCode;}
string Atom::get_charge() const {return charge;}
string Atom::get_resName() const {return resName;}
string Atom::get_name() const {return name;}
string Atom::get_element() const {return element;}
string Atom::get_recName() const {return "ATOM  ";}

double Atom::get_mass() const {
    if (__builtin_expect(element == "" || resName == "" || name == "", false)) {
        throw except::invalid_argument("Error in Atom::get_mass: Attempted to get atomic mass, but the element was not set!");
    }
    // mass of this nucleus + mass of attached H atoms
    return constants::mass::atomic.at(element) + constants::hydrogen_atoms::get.at(this->resName).at(this->name)*constants::mass::atomic.at("H");
};

void Atom::print() const {
    cout << "\nAtom no: " << serial << endl;
    cout << setw(17) << "(x, y, z): (" << setw(6) << coords.x() << ", " << setw(6) << coords.y() << ", " << setw(6) << coords.z() << ")" << endl;
    cout << setw(16) << "Weight: " << std::to_string(occupancy) << endl;
    cout << setw(16) << "Symbol: " << element << endl;
    cout << setw(16) << "Molecule: " << name << endl;
    return;
}

bool Atom::operator<(const Atom& rhs) const {
    return serial < rhs.serial;
}

bool Atom::operator==(const Atom& rhs) const {
    return uid == rhs.uid;
}

bool Atom::equals(const Atom& rhs) const {
    return operator==(rhs);
}

bool Atom::equals_content(const Atom& rhs) const {
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
    if (effective_charge != rhs.effective_charge) {return false;}
    return true;
}

Atom& Atom::operator=(const Atom& rhs) {
    name = rhs.name; 
    altLoc = rhs.altLoc; 
    resName = rhs.resName; 
    chainID = rhs.chainID; 
    iCode = rhs.iCode; 
    element = rhs.element; 
    charge = rhs.charge;
    occupancy = rhs.occupancy; 
    tempFactor = rhs.tempFactor;
    serial = rhs.serial; 
    resSeq = rhs.resSeq;
    coords = rhs.coords;
    effective_charge = rhs.effective_charge;
    uid = rhs.uid;
    return *this;
}