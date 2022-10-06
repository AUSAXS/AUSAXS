#pragma once

#include <string>
#include <vector>
#include <memory>

#include <data/Record.h>
#include <math/Vector3.h>

class Atom : public Record {
  public:
    /** 
     * @brief Construct a new Atom object.
     * 
     * @param v a TVector3 containing the x, y, z coordinates of the atom. 
     * @param occupancy the occupancy of this atom.
     * @param element the atomic element of the base atom.
     * @param name the molecule (e.g. HOH).
     * @param serial the serial number of this atom.
     */
    Atom(Vector3<double> v, double occupancy, std::string element, std::string name, int serial);

    /**
     * @brief Construct a new Atom object.
     * 
     * @param all see http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
     */
    Atom(int serial, std::string name, std::string altLoc, std::string resName, std::string chainID, int resSeq, std::string iCode, 
        Vector3<double> coords, double occupancy, double tempFactor, std::string element, std::string charge);

    /**
     * @brief Construct a new empty Atom object.
     */
    Atom();

    virtual ~Atom() override {}

    RecordType get_type() const override;

    /**
     * @brief Set the properties of this Atom based on a .pdb format string. 
     * 
     * @param s A .pdb format ATOM string. 
     */
    void parse_pdb(std::string s) override;

    /**
     * @brief Create a .pdb format string representation of this Atom. 
     * @return The string representation of this Atom. 
     */
    std::string as_pdb() const override;

    /** 
     * @brief Calculate the distance to another atom. 
     */
    double distance(const Atom& a) const;

    /** 
     * @brief Prints the contents of this object to the terminal. (NOT FULLY IMPLEMENTED!)
     */
    void print() const;

    /** 
     * @brief Move this atom by a vector.
     * 
     * @param v the translation vector.
     */
    void translate(Vector3<double> v);

    /**
     * @brief Determine if this is a water molecule.
     */
    virtual bool is_water() const;

//*** setters ***//
    void set_coordinates(Vector3<double> v);
    void set_x(double x);
    void set_y(double y);
    void set_z(double z);
    void set_occupancy(double occupancy);
    void set_tempFactor(double tempFactor);
    void set_altLoc(std::string altLoc);
    void set_serial(int serial);
    void set_resSeq(int resSeq);
    void set_effective_charge(double charge);
    void set_chainID(std::string chainID);
    void set_iCode(std::string iCode);
    void set_charge(std::string charge);

    /**
     * @brief Set the residue name for this atom.
     * 
     * @param resName the residue name, typically an amino acid such as LYS.
     */
    void set_resName(std::string resName);

    /**
     * @brief Specify the position of this atom within its residue.
     * 
     * @param name the position specifier, e.g. CG2 (Carbon | position G | branch 2).
     */
    void set_name(std::string name);

    /**
     * @brief Set the atomic element for this atom. Any spaces are removed. 
     * 
     * @param element the atomic element, e.g. He.
     */
    void set_element(std::string element);

//*** getters ***//
    Vector3<double>& get_coordinates();
    const Vector3<double>& get_coordinates() const;
    int get_serial() const;
    int get_resSeq() const;
    double get_occupancy() const;
    double get_tempFactor() const;
    double get_effective_charge() const;
    std::string get_altLoc() const;
    std::string get_chainID() const;
    std::string get_iCode() const;
    std::string get_charge() const;
    std::string get_resName() const;
    std::string get_name() const;
    std::string get_element() const;
    virtual std::string get_recName() const;

    /**
     * @brief Get the mass of this Atom.
     * 
     * @return The mass in u.
     */
    double get_mass() const;

    /**
     * @brief Get the number of protons in this atom.
     */
    unsigned int Z() const;

    /**
     * @brief Add @p charge to the effective charge of this atom. 
     */
    void add_effective_charge(const double charge) {effective_charge += charge;}

    /**
     * @brief Comparison function to allow this class to be a map key. 
     * 
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

    // properties as defined in https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, page 180.
    Vector3<double> coords = {0, 0, 0};
    std::string name, altLoc, resName, chainID, iCode, element, charge, recName;
    double occupancy = -1, tempFactor = -1;
    int serial = -1, resSeq = -1; 

    // other properties
    double effective_charge = -1;
    int uid = -1;

    // global counter for unique ids
    static inline int uid_counter = 0;
};