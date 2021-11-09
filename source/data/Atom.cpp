#pragma once

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
#include "Tools.cpp"

using namespace ROOT;
using std::vector, std::string, std::cout, std::endl, std::setw, std::left, std::right, std::shared_ptr, std::unique_ptr;
using boost::format;

class Atom : public Record {
public:
    /** 
     * @brief Constructor for the Atom class. 
     * @param v a TVector3 containing the x, y, z coordinates of the atom. 
     * @param occupancy the occupancy of this atom
     * @param element the atomic element of the base atom
     * @param name the molecule (e.g. HOH)
     */
    Atom(const TVector3 v, const double occupancy, const string element, const string name, int serial) {
        // we use our setters so we can validate the input if necessary
        this->set_coordinates(v);
        this->set_occupancy(occupancy);
        this->set_element(element);
        this->set_name(name);
        this->set_serial(serial);
    }

    Atom() {};

    RecordType get_type() override {return ATOM;}

    /**
     * @brief Set the properties of this Atom based on a .pdb format string. 
     * @param s A .pdb format ATOM string. 
     */
    void parse_pdb(const string ds) override {
        // THIS FORMAT DOES *NOT* FOLLOW THE STANDARD!
        // Modificiations: swapped space before 
        // http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
        // https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf

        string s = "ATOM  19152  NH2AARG A 281 12AB 14.421  84.308  88.373  0.50 40.45          He2+"; // DEBUG STRING
                //    "ATOM  19167  CB  PRO C 272     170.982 182.717 180.717  1.00 96.10           C  "
        //                    RE SE    NA AL RN  CH RS iC    X  Y  Z  OC TF     EL CH
        const char form[] = "%6c%5c%2c%4c%1c%3c %1c%4c%1c%4c%7c%8c%8c%6c%6c%10c%2c%2c";
        char recName[7] = "      ", serial[6] = "     ", space1[3] = "  ", name[5] = "    ", altLoc[2] = " ", resName[4] = "   ", chainID[2] = " ", resSeq[5] = "    ", 
            iCode[2] = " ", space2[5] = "    ", x[9] = "        ", y[9] = "        ", z[9] = "        ", occupancy[7] = "      ", tempFactor[7] = "      ", 
            space3[11] = "          ", element[3] = "  ", charge[3] = "  ";
        sscanf(s.c_str(), form, recName, serial, space1, name, altLoc, resName, chainID, resSeq, iCode, space2, x, y, z, occupancy, tempFactor, space3, element, charge);
        if (!(string(recName) == "ATOM  " || string(recName) == "HETATM")) {
            print_err("Error in Atom::parse_pdb: input string is not ATOM or HETATM (" + string(recName) + ").");
            exit(1);
        }

        // set all of the properties
        this->set_serial(std::atoi(serial));
        this->set_name(name);
        this->altLoc = string(altLoc);
        this->resName = string(resName);
        this->chainID = string(chainID);
        this->resSeq = std::atoi(resSeq);
        this->iCode = string(iCode);
        this->set_coordinates({std::atof(x), std::atof(y), std::atof(z)});
        this->set_occupancy(std::atof(occupancy));
        this->tempFactor = std::stod(tempFactor);
        this->set_element(element);
        this->charge = string(charge);

        // DEBUG OUTPUT
        cout << s << endl;
        cout << recName << serial << space1 << name << altLoc << resName << chainID << resSeq << iCode << space2 << x << y << z << occupancy << tempFactor << space3 << element << charge << endl;
        cout << this->as_pdb() << endl;
        cout << "name:" << name << ", altLoc:" << altLoc << ", resName:" << resName << ", chainID:" << chainID << ", resSeq:" << resSeq << ", iCode:" << iCode << endl;
        cout << "x:" << x << ", y:" << y << ", z" << z << ", occ:" << occupancy << ", tf:" << tempFactor << ", element:" << element << ", charge:" << charge << endl;
    }

    /**
     * @brief Create a .pdb format string representation of this Atom. 
     * @return The string representation of this Atom. 
     */
    string as_pdb() const override {
        std::stringstream ss;
        //                    RE SE    NA AL RN  CH RS iC    X  Y  Z  OC TF     EL CH
        const char form[] = "%6c%5c%2c%4c%1c%3c %1c%4c%1c%3c%8c%8c%8c%6c%6c%10c%2c%2c";
        ss << left << setw(6) << "ATOM  "                       // starts at index 0
            << right << setw(5) << get_serial() 
            << "  "          // 7
            << left << setw(4) << get_name()             // 13
            << left << setw(1) << altLoc                 // 17
            << right << setw(3) << resName
            << left << setw(1) << chainID                  // 23
            << " "
            << left << setw(4) << resSeq                                           // 27
            << right << setw(1) << iCode 
            << "  "
            << right << setw(8) << get_x()               // 31
            << right << setw(8) << get_y()               // 39
            << right << setw(8) << get_z()               // 47
            << right << setw(6) << get_occupancy()       // 55
            << right << setw(6) << tempFactor               // 61
            << "          "                                       // 67
            << left << setw(2) << get_element()          // 75
            << right << setw(2) << charge
            << endl;                                            // 79
        return ss.str();
    }

    /** Calculate the distance to another atom. 
     * @param a the other atom.
     * @return the distance. 
     */
    double distance(const shared_ptr<Atom> a) {
        return sqrt(pow(get_x() - a->get_x(), 2) + pow(get_y() - a->get_y(), 2) + pow(get_z() - a->get_z(), 2));
    }

    /** Prints the contents of this object to the terminal. */
    void print() {
        cout << "\nAtom no: " << serial << endl;
        cout << setw(17) << "(x, y, z): (" << setw(6) << get_x() << ", " << setw(6) << get_y() << ", " << setw(6) << get_z() << ")" << endl;
        cout << setw(16) << "Weight: " << std::to_string(occupancy) << endl;
        cout << setw(16) << "Symbol: " << element << endl;
        cout << setw(16) << "Molecule: " << name << endl;
        return;
    }

    /** Move this atom by a vector.
     * @param v the translation vector.
     */
    void translate(const TVector3 v) {
        coords += v;
    }

    /**
     * @brief Determine if this is a water molecule. 
     * @return true if this is a water molecule, otherwise false. 
     */
    bool is_water() {
        if (name == "HOH") {
            return true;
        }
        return false;
    }

    // setters
    void set_coordinates(TVector3 v) {this->coords = v;}
    void set_x(double x) {this->coords.SetX(x);}
    void set_y(double y) {this->coords.SetY(y);}
    void set_z(double z) {this->coords.SetZ(z);}
    void set_occupancy(double occupancy) {this->occupancy = occupancy;}
    void set_serial(int serial) {this->serial = serial;}

    /**
     * @brief Set the atomic name for this Atom. Any spaces are removed.
     * @param name the atomic name, e.g. LYS.
     */
    void set_name(string name) {
        boost::erase_all(name, " ");
        this->name = name;
    }

    /**
     * @brief Set the atomic element for this Atom. Any spaces are removed. 
     * @param element the atomic element.
     */
    void set_element(string element) {
        boost::erase_all(element, " ");
        if (!atomic_weight_map.count(element)) {
            print_err((format("ERROR: Invalid element \"%1%\".") % element).str());
            exit(1);
        }
        this->element = element;
    }

    // getters
    double get_x() const {return coords.X();}
    double get_y() const {return coords.Y();}
    double get_z() const {return coords.Z();}
    TVector3 get_coords() const {return coords;}
    double get_occupancy() const {return occupancy;}
    int get_serial() const {return serial;}
    string get_element() const {return element;}
    string get_name() const {return name;}

    double get_atomic_weight() const {
        if (element == "") {
            print_err("ERROR: Attempted to get atomic weight, but the element was not set!");
            exit(1);
        }
        return atomic_weight_map.at(element);
    };

    /**
     * @brief Comparison function to allow this class to be a map key. 
     * @param rhs Atom to compare against.
     */
    bool operator<(const Atom& rhs) const {
        return serial < rhs.get_serial();
    }

private:
    // properties as defined in https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, page 180.
    string name, altLoc, resName, chainID, iCode, element, charge;
    double occupancy, tempFactor;
    int serial, resSeq;
    TVector3 coords;

    // atomic weights taken from https://www.britannica.com/science/atomic-weight
    const std::map<string, int> atomic_weight_map = {{"H", 1.01}, {"He", 4.00}, {"Li", 6.95}, {"C", 12.01}, {"N", 14.01}, {"O", 16.00}, {"S", 32.06}};
};