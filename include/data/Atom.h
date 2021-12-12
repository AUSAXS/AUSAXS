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
#include "constants.h"
#include "math/Vector3.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using boost::format;

class Atom : public Record {
public:
    Atom(const Atom& a) : _name(a.name), _altLoc(a.altLoc), _resName(a.resName), _chainID(a.chainID), _iCode(a.iCode), _element(a.element), 
        _charge(a.charge), _occupancy(a.occupancy), _tempFactor(a.tempFactor), _serial(a.serial), _resSeq(a.resSeq), _coords(a.coords),
        _effective_charge(a.effective_charge), _uid(a.uid) {}

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
        set_coordinates(v);
        set_occupancy(occupancy);
        set_element(element);
        set_name(name);
        set_serial(serial);
        _altLoc = "";
        _resName = "";
        _chainID = "";
        _iCode = "";
        _charge = "";
        _resSeq = -1;
        _tempFactor = -1;
        _uid = uid_counter++;
    }

    /**
     * @brief Construct a new Atom object.
     * @param all see http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
     */
    Atom(const int serial, const string name, const string altLoc, const string resName, const string chainID, const int resSeq, 
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
            _effective_charge = constants::charge::get.at(this->element) + constants::hydrogen_atoms::get.at(this->resName).at(this->name);
            _uid = uid_counter++;
    }

    /**
     * @brief Construct a new empty Atom object.
     */
    Atom() {
        _name = "";
        _altLoc = "";
        _resName = "";
        _chainID = "";
        _iCode = "";
        _element = "";
        _charge = "";
        _serial = -1;
        _resSeq = -1;
        _occupancy = -1;
        _tempFactor = -1;
        _coords = {0, 0, 0};
        _effective_charge = -1;
        _uid = uid_counter++;
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

    /** 
     * @brief Calculate the distance to another atom. 
     */
    double distance(const Atom& a) const {return coords.distance(a.coords);}

    /** 
     * @brief Prints the contents of this object to the terminal. (NOT FULLY IMPLEMENTED!)
     */
    void print();

    /** 
     * @brief Move this atom by a vector.
     * @param v the translation vector.
     */
    void translate(const Vector3 v) {
        _coords += v;
    }

    /**
     * @brief Determine if this is a water molecule. Only used by the Hetatom subclass, but is defined here for convenience. 
     * @return true if this is a water molecule, otherwise false. 
     */
    virtual bool is_water() const {return false;}

    // setters
    void set_coordinates(const Vector3 v) {_coords = v;}
    void set_x(const double x) {_coords.x = x;}
    void set_y(const double y) {_coords.y = y;}
    void set_z(const double z) {_coords.z = z;}
    void set_occupancy(double occupancy) {_occupancy = occupancy;}
    void set_serial(const int serial) {_serial = serial;}
    void set_resSeq(const int resSeq) {_resSeq = resSeq;}
    void set_effective_charge(double charge) {_effective_charge = charge;}

    /**
     * @brief Set the residue name for this atom.
     * @param resName the residue name, typically an amino acid such as LYS.
     */
    void set_resName(string resName) {_resName = resName;}
    void set_chainID(string chainID) {_chainID = chainID;}

    /**
     * @brief Specify the position of this atom within its residue.
     * @param name the position specifier, e.g. CG2 (Carbon | position G | branch 2).
     */
    void set_name(string name) {_name = name;}

    /**
     * @brief Set the atomic element for this atom. Any spaces are removed. 
     * @param element the atomic element, e.g. He.
     */
    void set_element(string element) {
        if (__builtin_expect(constants::mass::atomic.count(element) == 0, false)) { // check that the weight is defined
            print_err((format("Error in Atom::set_element: The weight of element \"%1%\" is not defined.") % element).str());
            exit(1);
        }
        _element = element;
    }

    virtual string get_recName() const {return "ATOM  ";}

    double get_mass() const {
        if (__builtin_expect(element == "" || resName == "" || name == "", false)) {
            print_err("Error in Atom::get_atomic_weight: Attempted to get atomic mass, but the element was not set!");
            exit(1);
        }
        // mass of this nucleus + mass of attached H atoms
        return constants::mass::atomic.at(element) + constants::hydrogen_atoms::get.at(this->resName).at(this->name)*constants::mass::atomic.at("H");
    };

    /**
     * @brief Add @p charge to the effective charge of this atom. 
     */
    void add_effective_charge(const double& charge) {_effective_charge += charge;}

    /**
     * @brief Comparison function to allow this class to be a map key. 
     * @param rhs Atom to compare against.
     */
    bool operator<(const Atom& rhs) const;

    /**
     * @brief Equality operator to determine if two atoms are equal
     *        Note that this is a @a content comparator, and thus determines if two atoms are equal based on their contents. 
     * @param rhs the atom to compare against. 
     */
    bool operator==(const Atom& rhs) const;

    Atom& operator=(const Atom& rhs) {
        _name = rhs.name; _altLoc = rhs.altLoc; _resName = rhs.resName; _chainID = rhs.chainID; _iCode = rhs.iCode; _element = rhs.element; _charge = rhs.charge;
        _occupancy = rhs.occupancy; _tempFactor = rhs.tempFactor;
        _serial = rhs.serial; _resSeq = rhs.resSeq;
        _coords = rhs.coords;
        _effective_charge = rhs.effective_charge;
        _uid = rhs.uid;
        return *this;
    }

    const string &name = _name, &altLoc = _altLoc, &resName = _resName, &chainID = _chainID, &iCode = _iCode, &element = _element, &charge = _charge;
    const double &occupancy = _occupancy, &tempFactor = _tempFactor;
    const int &serial = _serial, &resSeq = _resSeq;
    const Vector3& coords = _coords;
    const double& effective_charge = _effective_charge;
    const int &uid = _uid;
protected:
    // properties as defined in https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, page 180.
    string _name, _altLoc, _resName, _chainID, _iCode, _element, _charge;
    double _occupancy, _tempFactor;
    int _serial, _resSeq; 
    Vector3 _coords;

    // other properties
    double _effective_charge;

    int _uid;
    static inline int uid_counter = 0;
};