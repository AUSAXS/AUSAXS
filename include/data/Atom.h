#pragma once

#include <map>
#include <string>
#include <vector>
#include <utility>

#include "data/Record.h"
#include "Tools.h"
#include "constants.h"
#include "math/Vector3.h"
#include "Exceptions.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;

class Atom : public Record {
  public:
    Atom(const Atom&& a) noexcept;

    Atom(const Atom& a);

    /** 
     * @brief Construct a new Atom object.
     * @param v a TVector3 containing the x, y, z coordinates of the atom. 
     * @param occupancy the occupancy of this atom.
     * @param element the atomic element of the base atom.
     * @param name the molecule (e.g. HOH).
     * @param serial the serial number of this atom.
     */
    Atom(const Vector3 v, const double occupancy, const string element, const string name, int serial);

    /**
     * @brief Construct a new Atom object.
     * @param all see http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
     */
    Atom(const int serial, const string name, const string altLoc, const string resName, const string chainID, const int resSeq, 
        const string iCode, const Vector3 coords, const double occupancy, const double tempFactor, const string element, const string charge);

    /**
     * @brief Construct a new empty Atom object.
     */
    Atom();

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
    void print() const;

    /** 
     * @brief Move this atom by a vector.
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



//*** setters ***//
    void set_coordinates(const Vector3 v) {coords = v;}
    void set_x(double x) {coords.x = x;}
    void set_y(double y) {coords.y = y;}
    void set_z(double z) {coords.z = z;}
    void set_occupancy(double occupancy) {this->occupancy = occupancy;}
    void set_tempFactor(double tempFactor) {this->tempFactor = tempFactor;}
    void set_altLoc(string altLoc) {this->altLoc = altLoc;}
    void set_serial(int serial) {this->serial = serial;}
    void set_resSeq(int resSeq) {this->resSeq = resSeq;}
    void set_effective_charge(double charge) {effective_charge = charge;}
    void set_chainID(string chainID) {this->chainID = chainID;}
    void set_iCode(string iCode) {this->iCode = iCode;}
    void set_charge(string charge) {this->charge = charge;}

    /**
     * @brief Set the residue name for this atom.
     * @param resName the residue name, typically an amino acid such as LYS.
     */
    void set_resName(string resName) {this->resName = resName;}

    /**
     * @brief Specify the position of this atom within its residue.
     * @param name the position specifier, e.g. CG2 (Carbon | position G | branch 2).
     */
    void set_name(string name) {this->name = name;}

    /**
     * @brief Set the atomic element for this atom. Any spaces are removed. 
     * @param element the atomic element, e.g. He.
     */
    void set_element(string element) {
        if (__builtin_expect(constants::mass::atomic.count(element) == 0, false)) { // check that the weight is defined
            throw except::invalid_argument("Error in Atom::set_element: The weight of element " + element + " is not defined.");
        }
        this->element = element;
    }

//*** getters ***//
    Vector3& get_coordinates() {return coords;}
    const Vector3& get_coordinates() const {return coords;}
    int get_serial() const {return serial;}
    int get_resSeq() const {return resSeq;}
    double get_occupancy() const {return occupancy;}
    double get_tempFactor() const {return tempFactor;}
    double get_effective_charge() const {return effective_charge;}
    string get_altLoc() const {return altLoc;}
    string get_chainID() const {return chainID;}
    string get_iCode() const {return iCode;}
    string get_charge() const {return charge;}
    string get_resName() const {return resName;}
    string get_name() const {return name;}
    string get_element() const {return element;}
    virtual string get_recName() const {return "ATOM  ";}

    double get_mass() const {
        if (__builtin_expect(element == "" || resName == "" || name == "", false)) {
            throw except::invalid_argument("Error in Atom::get_atomic_weight: Attempted to get atomic mass, but the element was not set!");
        }
        // mass of this nucleus + mass of attached H atoms
        return constants::mass::atomic.at(element) + constants::hydrogen_atoms::get.at(this->resName).at(this->name)*constants::mass::atomic.at("H");
    };

    /**
     * @brief Add @p charge to the effective charge of this atom. 
     */
    void add_effective_charge(const double charge) {effective_charge += charge;}

    /**
     * @brief Comparison function to allow this class to be a map key. 
     * @param rhs Atom to compare against.
     */
    bool operator<(const Atom& rhs) const;

    /**
     * @brief Equality operator to determine if two atoms are equal.
     *        Note that this compares their unique object identifier which is generated at object creation, completely disregarding
     *        their contents. Unless a deliberate attempt at desyncing the id from the contents were made, equality of content follows
     *        from equality of id. 
     * @param rhs Atom to compare against. 
     */
    bool operator==(const Atom& rhs) const;

    /**
     * @brief Equality operator to determine if two atoms are equal.
     *        Note that this is a @a content comparator, and thus determines if two atoms are equal based on their contents. 
     * @param rhs Atom to compare against. 
     */
    bool equals_content(const Atom& rhs) const;

    /**
     * @brief Equality operator to determine if two atoms are equal.
     *        Note that this compares their unique object identifier which is generated at object creation, completely disregarding
     *        their contents. Unless a deliberate attempt at desyncing the id from the contents were made, equality of content follows
     *        from equality of id. 
     * @param rhs Atom to compare against. 
     */
    bool equals(const Atom& rhs) const;

    /**
     * @brief Inequality operator to determine if two atoms are not equal.
     *        Note that this compares their unique object identifier which is generated at object creation, completely disregarding
     *        their contents. Unless a deliberate attempt at desyncing the id from the contents were made, inequality of content follows
     *        from inequality of id. 
     * @param rhs Atom to compare against. 
     */
    bool operator!=(const Atom& rhs) const {return !operator==(rhs);}

    Atom& operator=(const Atom& rhs);

    // properties as defined in https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, page 180.
    Vector3 coords = {0, 0, 0};
    string name = "", altLoc = "", resName = "", chainID = "", iCode = "", element = "", charge = "";
    double occupancy = -1, tempFactor = -1;
    int serial = -1, resSeq = -1; 

    // other properties
    double effective_charge = -1;
    int uid = -1;

    // global counter for unique ids
    static inline int uid_counter = 0;
};