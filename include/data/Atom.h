#pragma once

#include <data/Record.h>
#include <math/Vector3.h>

#include <string>
#include <vector>
#include <memory>

namespace constants {enum class atom_t;}
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
        Atom(Vector3<double> v, double occupancy, constants::atom_t element, const std::string& name, int serial);

        /**
         * @brief Construct a new Atom object.
         * 
         * @param all see http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
         */
        Atom(int serial, const std::string& name, const std::string& altLoc, const std::string& resName, char chainID, int resSeq, 
            const std::string& iCode, Vector3<double> coords, double occupancy, double tempFactor, constants::atom_t element, const std::string& charge);

        /**
         * @brief Default constructor.
         */
        Atom();

        /**
         * @brief Destructor.
         */
        virtual ~Atom() override {}

        RecordType get_type() const override;

        /**
         * @brief Set the properties of this Atom based on a .pdb format ATOM string. 
         */
        void parse_pdb(const std::string& s) override;

        /**
         * @brief Create a .pdb format string representation of this Atom. 
         */
        std::string as_pdb() const override;

        /** 
         * @brief Calculate the distance to another atom. 
         */
        double distance(const Atom& a) const;

        /** 
         * @brief Translate this atom.
         */
        void translate(Vector3<double> v);

        /**
         * @brief Check if this is a water molecule by checking the residue name. Both HOH and SOL are considered water molecules.
         */
        virtual bool is_water() const;

    //*** setters ***//
        /**
         * @brief Set the atomic coordinates.
         */
        void set_coordinates(Vector3<double> v);


        /**
         * @brief Set the x coordinate.
         */
        void set_x(double x);

        /**
         * @brief Set the y coordinate.
         */
        void set_y(double y);

        /**
         * @brief Set the z coordinate.
         */
        void set_z(double z);

        /**
         * @brief Set the occupancy.
         * 
         * @param occupancy the occupancy.
         */
        void set_occupancy(double occupancy);

        /**
         * @brief Set the temperature factor for this atom.
         *        Currently this is not used for anything. 
         * 
         * @param tempFactor the temperature factor.
         */
        void set_tempFactor(double tempFactor);

        /**
         * @brief Set the alternate location for this atom.
         * 
         * @param altLoc the alternate location, e.g. A.
         */
        void set_altLoc(const std::string& altLoc);

        /**
         * @brief Set the serial number for this atom.
         * 
         * @param serial the serial number, e.g. 1.
         */
        void set_serial(int serial);

        /**
         * @brief Set the residue sequence number for this atom.
         * 
         * @param resSeq the residue sequence number, e.g. 1.
         */
        void set_resSeq(int resSeq);

        /**
         * @brief Set the effective charge for this atom.
         * 
         * @param charge the effective charge, e.g. +1.
         */
        void set_effective_charge(double charge);

        /**
         * @brief Set the chain ID for this atom.
         * 
         * @param chainID the chain ID, e.g. A.
         */
        void set_chainID(char chainID);

        /**
         * @brief Set the insertion code for this atom.
         * 
         * @param iCode the insertion code, e.g. A.
         */
        void set_iCode(const std::string& iCode);

        /**
         * @brief Set the charge for this atom.
         * 
         * @param charge the charge, e.g. +1.
         */
        void set_charge(const std::string& charge);

        /**
         * @brief Set the residue name for this atom.
         * 
         * @param resName the residue name, typically an amino acid such as LYS.
         */
        void set_resName(const std::string& resName);

        /**
         * @brief Specify the position of this atom within its residue.
         * 
         * @param name the position specifier, e.g. CG2 (Carbon | position G | branch 2).
         */
        void set_name(const std::string& name);

        /**
         * @brief Set the atomic element for this atom. Any spaces are removed. 
         * 
         * @param element the atomic element, e.g. He.
         */
        void set_element(constants::atom_t element);

        /**
         * @brief Set the atomic element for this atom. Any spaces are removed. 
         * 
         * @param element the atomic element, e.g. He.
         */
        void set_element(const std::string& element);

        /**
         * @brief Get the coordinates of this Atom.
         */
        Vector3<double>& get_coordinates();

        /**
         * @brief Get the coordinates of this Atom.
         */
        const Vector3<double>& get_coordinates() const;

        /**
         * @brief Get the serial of this Atom.
         */
        int get_serial() const;

        /**
         * @brief Get the residue sequence number of this Atom.
         */
        int get_resSeq() const;

        /**
         * @brief Get the occupancy of this Atom.
         */
        double get_occupancy() const;

        /**
         * @brief Get the temperature factor of this Atom.
         */
        double get_tempFactor() const;

        /**
         * @brief Get the effective charge of this Atom. 
         *        The effective charge is the sum of the absolute charge and the number of bound hydrogens.
         *        Note that this charge can be modified during runtime.
         */
        double get_effective_charge() const;

        /**
         * @brief Get the absolute charge of this Atom. The absolute charge is just the Z value of the atom. (Equivalent to Atom::Z)
         */
        double get_absolute_charge() const;

        /**
         * @brief Get the alternate location of this Atom.
         */
        std::string get_altLoc() const;

        /**
         * @brief Get the chain ID of this Atom.
         */
        char get_chainID() const;

        /**
         * @brief Get the insertion code of this Atom.
         */
        std::string get_iCode() const;

        /**
         * @brief Get the string representation of the charge of this Atom.
         */
        std::string get_charge() const;

        /**
         * @brief Get the residue name of this Atom.
         */
        std::string get_resName() const;

        /**
         * @brief Get the name of this Atom.
         */
        std::string get_name() const;

        /**
         * @brief Get the atomic element of this Atom.
         */
        constants::atom_t get_element() const;

        /**
         * @brief Get the record name of this Atom. 
         *        This is the first 6 characters of the line in the PDB file.
         */
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
        void add_effective_charge(const double charge);

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
         * @brief Inequality operator to determine if two atoms are not equal.
         *        Note that this compares their unique object identifier which is generated at object creation, completely disregarding
         *        their contents. Unless a deliberate attempt at desyncing the id from the contents were made, inequality of content follows
         *        from inequality of id. 
         * @param rhs Atom to compare against. 
         */
        bool operator!=(const Atom& rhs) const {return !operator==(rhs);}

        // properties as defined in https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, page 180.
        Vector3<double> coords = {0, 0, 0};
        std::string name, altLoc, resName, iCode, charge, recName = "ATOM  ";
        char chainID;
        constants::atom_t element;
        double occupancy = -1, tempFactor = -1;
        int serial = -1, resSeq = -1; 

        // other properties
        double effective_charge = -1;
        int uid = -1;

    private: 
        static inline int uid_counter = 0; // global counter for unique ids
};