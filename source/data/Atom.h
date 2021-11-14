#pragma once

// includes
#include <map>
#include <string>
#include <vector>
#include <utility>
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
    void parse_pdb(const string s) override;

    /**
     * @brief Create a .pdb format string representation of this Atom. 
     * @return The string representation of this Atom. 
     */
    string as_pdb() const override;

    /** Calculate the distance to another atom. 
     * @param a the other atom.
     * @return the distance. 
     */
    double distance(const shared_ptr<Atom> a) {
        return sqrt(pow(get_x() - a->get_x(), 2) + pow(get_y() - a->get_y(), 2) + pow(get_z() - a->get_z(), 2));
    }

    /** Prints the contents of this object to the terminal. */
    void print();

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
    bool is_water() const;

    // setters
    void set_coordinates(const TVector3 v) {this->coords = v;}
    void set_x(const double x) {this->coords.SetX(x);}
    void set_y(const double y) {this->coords.SetY(y);}
    void set_z(const double z) {this->coords.SetZ(z);}
    void set_occupancy(double occupancy) {this->occupancy = occupancy;}
    void set_serial(const int serial) {this->serial = serial;}
    void set_resSeq(const int resSeq) {this->resSeq = resSeq;}
    void set_resName(string resName) {
        boost::erase_all(resName, " "); // remove spaces
        this->resName = resName;
    }
    void set_chainID(string chainID) {
        boost::erase_all(chainID, " "); // remove spaces
        this->chainID = chainID;
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
    int get_resSeq() const {return resSeq;}
    int get_serial() const {return serial;}
    double get_occupancy() const {return occupancy;}
    double get_x() const {return coords.X();}
    double get_y() const {return coords.Y();}
    double get_z() const {return coords.Z();}
    TVector3 get_coords() const {return coords;}
    string get_element() const {return element;}
    string get_resName() const {return resName;}
    string get_iCode() const {return iCode;}
    string get_chainID() const {return chainID;}

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
    bool operator<(const Atom& rhs) const;

    /**
     * @brief Equality operator to determine if two atoms are equal
     * @param rhs the atom to compare against. 
     * @return True if equal, false otherwise.
     */
    bool operator==(const Atom& rhs) const;

    /**
     * @brief Create a new water Atom.
     * @return A pointer to the new water Atom. 
     */
    static unique_ptr<Atom> create_new_water();

    /**
     * @brief Create a new water Atom.
     * @param coords the coordinates for the new Atom.
     * @return A pointer to the new water Atom. 
     */
    static unique_ptr<Atom> create_new_water(TVector3 coords);

private:
    // properties as defined in https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, page 180.
    string recName, name, altLoc, resName, chainID, iCode, element, charge;
    double occupancy, tempFactor;
    int serial, resSeq; 
    TVector3 coords;

    // atomic weights taken from https://www.britannica.com/science/atomic-weight
    const std::map<string, int> atomic_weight_map = {{"H", 1.01}, {"He", 4.00}, {"Li", 6.95}, {"C", 12.01}, {"N", 14.01}, {"O", 16.00}, {"S", 32.06}};
};