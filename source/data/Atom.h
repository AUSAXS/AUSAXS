#pragma once

// includes
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <boost/algorithm/string.hpp>

// my own stuff
#include "data/Record.h"
#include "Tools.h"
#include "data/properties.h"
#include "math/Vector3.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
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
    Atom(const Vector3 v, const double occupancy, const string element, const string name, int serial) {
        // we use our setters so we can validate the input if necessary
        this->set_coordinates(v);
        this->set_occupancy(occupancy);
        this->set_element(element);
        this->set_name(name);
        this->set_serial(serial);
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
    Atom(const int serial, const string name, const string altLoc, const string resName, const string chainID, const int resSeq, 
        const string iCode, const Vector3 coords, const double occupancy, const double tempFactor, const string element, const string charge) {
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

    virtual ~Atom() override {}

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
    void translate(const Vector3 v) {
        coords += v;
    }

    /**
     * @brief Determine if this is a water molecule. Only used by the Hetatom subclass, but is defined here for convenience. 
     * @return true if this is a water molecule, otherwise false. 
     */
    virtual bool is_water() const {return false;}

    // setters
    void set_coordinates(const Vector3 v) {this->coords = v;}
    void set_x(const double x) {this->coords.x = x;}
    void set_y(const double y) {this->coords.y = y;}
    void set_z(const double z) {this->coords.z = z;}
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
        if (!property::weight::atomic.count(element)) { // check that the weight is defined
            print_err((format("Error in Atom::set_element: The weight of element \"%1%\" is not defined.") % element).str());
            exit(1);
        }
        this->element = element;
    }

    // getters
    int get_resSeq() const {return resSeq;}
    int get_serial() const {return serial;}
    double get_occupancy() const {return occupancy;}
    double get_x() const {return coords.x;}
    double get_y() const {return coords.y;}
    double get_z() const {return coords.z;}
    Vector3 get_coords() const {return coords;}
    string get_element() const {return element;}
    string get_resName() const {return resName;}
    string get_iCode() const {return iCode;}
    string get_chainID() const {return chainID;}
    virtual string get_recName() const {return "ATOM  ";}

    double get_atomic_weight() const {
        if (element == "") {
            print_err("Error in Atom::get_atomic_weight: Attempted to get atomic weight, but the element was not set!");
            exit(1);
        }
        return property::weight::atomic.at(element);
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

protected:
    // properties as defined in https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, page 180.
    string name, altLoc, resName, chainID, iCode, element, charge;
    double occupancy, tempFactor;
    int serial, resSeq; 
    Vector3 coords;

    // other properties
    double effective_charge;
};