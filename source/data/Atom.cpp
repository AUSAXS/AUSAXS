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
     * @brief Construct a new Atom object.
     * @param v a TVector3 containing the x, y, z coordinates of the atom. 
     * @param occupancy the occupancy of this atom.
     * @param element the atomic element of the base atom.
     * @param name the molecule (e.g. HOH).
     * @param serial the serial number of this atom.
     */
    Atom(const TVector3 v, const double occupancy, const string element, const string name, int serial) {
        // we use our setters so we can validate the input if necessary
        this->set_coordinates(v);
        this->set_occupancy(occupancy);
        this->set_element(element);
        this->set_name(name);
        this->set_serial(serial);
        this->recName = "";
        this->altLoc = "";
        this->resName = "";
        this->chainID = "";
        this->iCode = "";
        this->charge = "";
        this->resSeq = -1;
        this->tempFactor = -1;
    }

    /**
     * @brief Construct a new Atom object.
     * @param all see http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
     */
    Atom(const string recName, const int serial, const string name, const string altLoc, const string resName, const string chainID, const int resSeq, 
        const string iCode, const TVector3 coords, const double occupancy, const double tempFactor, const string element, const string charge) {
            this->recName = recName;
            this->set_serial(serial);
            this->set_name(name);
            this->altLoc = altLoc;
            this->set_resName(resName);
            this->chainID = chainID;
            this->resSeq = resSeq;
            this->iCode = iCode;
            this->set_coordinates(coords);
            this->set_occupancy(occupancy);
            this->tempFactor = tempFactor;
            this->set_element(element);
            this->charge = charge;
    }

    /**
     * @brief Construct a new empty Atom object.
     */
    Atom() {
        recName = "NULL";
        name = "";
        altLoc = "";
        resName = "";
        chainID = "";
        iCode = "";
        element = "";
        charge = "";
        serial = -1;
        resSeq = -1;
        occupancy = -1;
        tempFactor = -1;
        coords = {0, 0, 0};
    }

    RecordType get_type() const override {return ATOM;}

    /**
     * @brief Set the properties of this Atom based on a .pdb format string. 
     * @param s A .pdb format ATOM string. 
     */
    void parse_pdb(const string s) override {
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

    /**
     * @brief Create a .pdb format string representation of this Atom. 
     * @return The string representation of this Atom. 
     */
    string as_pdb() const override {
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
            << right << setw(8) << std::setprecision(7) << get_x()               // 31 - 38
            << right << setw(8) << std::setprecision(7) << get_y()               // 39 - 46
            << right << setw(8) << std::setprecision(7) << get_z()               // 47 - 54
            << right << setw(6) << std::setprecision(7) << get_occupancy()       // 55 - 60
            << right << setw(6) << std::setprecision(7) << tempFactor            // 61 - 66
            << "          "                                                      // 67 - 76
            << right << setw(2) << get_element()                                 // 77 - 78
            << left << setw(2) << charge                                         // 79 - 80
            << endl;
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
    bool is_water() const {
        if (resName == "HOH") {
            return true;
        }
        return false;
    }

    // setters
    void set_coordinates(const TVector3 v) {this->coords = v;}
    void set_x(const double x) {this->coords.SetX(x);}
    void set_y(const double y) {this->coords.SetY(y);}
    void set_z(const double z) {this->coords.SetZ(z);}
    void set_occupancy(double occupancy) {this->occupancy = occupancy;}
    void set_serial(const int serial) {this->serial = serial;}
    void set_resName(string resName) {
        boost::erase_all(resName, " "); // remove spaces
        this->resName = resName;
    }

    /**
     * @brief Set the atomic name for this Atom. Any spaces are removed.
     * @param name the atomic name, e.g. LYS.
     */
    void set_name(string name) {
        boost::erase_all(name, " "); // remove spaces
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
    string get_resName() const {return resName;}

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

    bool operator==(const Atom& rhs) const {
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

    /**
     * @brief Create a new water Atom.
     * @return A pointer to the new water Atom. 
     */
    static unique_ptr<Atom> create_new_water() {
        return create_new_water({0, 0, 0});
    }

    /**
     * @brief Create a new water Atom.
     * @param coords the coordinates for the new Atom.
     * @return A pointer to the new water Atom. 
     */
    static unique_ptr<Atom> create_new_water(TVector3 coords) {
        return std::make_unique<Atom>(Atom("HETATM", -1, "O", "", "HOH", "", -1, "", coords, 1, 0, "O", ""));
    }

private:
    // properties as defined in https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, page 180.
    string recName, name, altLoc, resName, chainID, iCode, element, charge;
    double occupancy, tempFactor;
    int serial, resSeq;
    TVector3 coords;

    // atomic weights taken from https://www.britannica.com/science/atomic-weight
    const std::map<string, int> atomic_weight_map = {{"H", 1.01}, {"He", 4.00}, {"Li", 6.95}, {"C", 12.01}, {"N", 14.01}, {"O", 16.00}, {"S", 32.06}};
};