#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iomanip>
#include <iostream>

#include <data/Record.h>
#include <data/Atom.h>
#include <math/Vector3.h>
#include <utility/Constants.h>
#include <utility/Settings.h>
#include <utility/Utility.h>

using std::vector, std::string, std::cout, std::endl, std::setw, std::left, std::right, std::shared_ptr, std::unique_ptr;

Atom::Atom(Vector3<double> v, double occupancy, string element, string resName, int serial) : uid(uid_counter++) {
    // we use our setters so we can validate the input if necessary
    set_coordinates(v);
    set_occupancy(occupancy);
    set_element(element);
    set_resName(resName);
    set_serial(serial);
    set_effective_charge(constants::charge::atomic.get(this->element));
}

Atom::Atom(int serial, string name, string altLoc, string resName, string chainID, int resSeq, string iCode, 
    Vector3<double> coords, double occupancy, double tempFactor, string element, string charge) : uid(uid_counter++) {
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

        // check if we want to correct the hydrogen contribution to the charge
        if (setting::protein::use_effective_charge) {
            // use a try-catch block to throw more sensible errors
            try {
                effective_charge = constants::charge::atomic.get(this->element) + constants::hydrogen_atoms::residues.get(this->resName).get(this->name, this->element);
            } catch (const except::base& e) {
                throw except::invalid_argument("Atom::Atom: Could not set effective charge. Unknown element, residual or atom: (" + element + ", " + resName + ", " + name + ")");
            }
        } 
        // otherwise just use the atomic charge (this will never fail)
        else {
            effective_charge = constants::charge::atomic.get(this->element);
        }
        uid = uid_counter++;
}

Atom::Atom() : uid(uid_counter++) {}

void Atom::parse_pdb(string s) {
    s = utility::remove_all(s, "\n\r"); // remove any newline or carriage return
    int pad_size = 81 - s.size();
    if (pad_size < 0) {
        utility::print_warning("Warning in Atom::parse_pdb: Line is longer than 80 characters. Truncating.");
        std::cout << "\"" << s << "\"" << std::endl;
        s = s.substr(0, 80);
    } else {
        s += string(pad_size, ' ');
    }

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
    if (!(Record::get_type(recName) == RecordType::ATOM)) [[unlikely]] {
        throw except::parse_error("Atom::parse_pdb: input string is not \"ATOM  \" or \"HETATM\" (" + recName + ").");
    }

    // remove any spaces from the numbers
    serial = utility::remove_all(serial, " ");
    name = utility::remove_all(space1+name, " "); // we add space1 since some programs (gromacs) uses it for the name.
    resName = utility::remove_all(resName, " ");
    resSeq = utility::remove_all(resSeq, " ");
    x = utility::remove_all(x, " ");
    y = utility::remove_all(y, " ");
    z = utility::remove_all(z, " ");
    occupancy = utility::remove_all(occupancy, " ");
    tempFactor = utility::remove_all(tempFactor, " ");
    element = utility::remove_all(element, " ");

    // sometimes people use the first character of x for some other purpose.
    // if it is a digit, the following won't work. On the other hand they're kinda asking for it then. Follow the standard, people. 
    if (!(std::isdigit(x[0]) || x[0] == '-')) {
        x = x.substr(1, x.size()-1);
    }

    // set all of the properties
    try {
        this->recName = recName;
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
    } catch (const except::base& e) { // catch conversion errors and output a more meaningful error message
        utility::print_warning("Atom::parse_pdb: Invalid field values in line \"" + s + "\".");
        throw;
    }

    if (setting::protein::use_effective_charge) {
        effective_charge = constants::charge::atomic.get(this->element) + constants::hydrogen_atoms::residues.get(this->resName).get(this->name, this->element);
    } else {
        effective_charge = constants::charge::atomic.get(this->element);
    }
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
        << right << setw(8) << utility::fixedwidth(coords.x(), 7)            // 31 - 38
        << right << setw(8) << utility::fixedwidth(coords.y(), 7)            // 39 - 46
        << right << setw(8) << utility::fixedwidth(coords.z(), 7)            // 47 - 54
        << right << setw(6) << utility::fixedwidth(occupancy, 6)             // 55 - 60
        << right << setw(6) << utility::fixedwidth(tempFactor, 6)            // 61 - 66
        << "          "                                                      // 67 - 76
        << right << setw(2) << element                                       // 77 - 78
        << left << setw(2) << charge                                         // 79 - 80
        << endl;
    return ss.str();
}

Record::RecordType Atom::get_type() const {return RecordType::ATOM;}

double Atom::distance(const Atom& a) const {return coords.distance(a.coords);}
void Atom::translate(Vector3<double> v) {coords += v;}
bool Atom::is_water() const {return (resName == "HOH") | (resName == "SOL");}
void Atom::set_coordinates(Vector3<double> v) {coords = v;}
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
    if (!constants::mass::atomic.contains(element)) [[unlikely]] { // check that the weight is defined
        throw except::invalid_argument("Atom::set_element: The weight of element " + element + " is not defined.");
    }
    this->element = element;
}

Vector3<double>& Atom::get_coordinates() {return coords;}
const Vector3<double>& Atom::get_coordinates() const {return coords;}
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
string Atom::get_recName() const {return recName;}

double Atom::get_mass() const {
    if (setting::protein::use_effective_charge) {
        // mass of this nucleus + mass of attached H atoms
        if (element.empty() || resName.empty() || name.empty()) [[unlikely]] {
            throw except::invalid_argument("Atom::get_mass: Attempted to get atomic mass, but the element, residue name, or name was not set!");
        }
        return constants::mass::atomic.get(element) + constants::hydrogen_atoms::residues.get(this->resName).get(this->name, this->element);
    } else {
        if (element.empty()) [[unlikely]] {
            throw except::invalid_argument("Atom::get_mass: Attempted to get atomic mass, but the element was not set!");
        }
        return constants::mass::atomic.get(element);
    }
}

unsigned int Atom::Z() const {
    if (element == "") [[unlikely]] {
        throw except::invalid_argument("Atom::get_Z: Attempted to get atomic charge, but the element was not set!");
    }
    return constants::charge::atomic.get(element);
}

void Atom::add_effective_charge(const double charge) {effective_charge += charge;}

bool Atom::operator<(const Atom& rhs) const {
    return serial < rhs.serial;
}

bool Atom::operator==(const Atom& rhs) const {
    return uid == rhs.uid;
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