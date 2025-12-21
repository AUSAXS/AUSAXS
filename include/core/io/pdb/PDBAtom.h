// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/pdb/Record.h>
#include <math/Vector3.h>
#include <constants/ConstantsFwd.h>
#include <form_factor/FormFactorType.h>

#include <string>

namespace ausaxs::io::pdb {
    class PDBAtom : public Record {
        public:
			PDBAtom();
			PDBAtom(const PDBAtom& rhs) = default;
			PDBAtom(PDBAtom&& rhs) noexcept = default;
			PDBAtom &operator=(const PDBAtom& rhs) = default;
			PDBAtom &operator=(PDBAtom&& rhs) noexcept = default;
			~PDBAtom() override = default;

            /** 
             * @brief Construct a new Atom object.
             * 
             * @param v a TVector3 containing the x, y, z coordinates of the atom. 
             * @param occupancy the occupancy of this atom.
             * @param element the atomic element of the base atom.
             * @param name the molecule (e.g. HOH).
             * @param serial the serial number of this atom.
             */
            PDBAtom(Vector3<double> v, double occupancy, constants::atom_t element, const std::string& name, int serial);

            /**
             * @brief Construct a new Atom object.
             * 
             * @param all see http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
             */
            PDBAtom(int serial, const std::string& name, const std::string& altLoc, const std::string& resName, char chainID, int resSeq, 
                const std::string& iCode, Vector3<double> coords, double occupancy, double tempFactor, constants::atom_t element, const std::string& charge);

            RecordType get_type() const override;

            form_factor::form_factor_t get_form_factor_type() const;

            /**
             * @brief Add implicit hydrogens to this atom. 
             *        This is done by adding the number of implicit hydrogens to the effective charge of the atom and modifying its form factor. 
             */
            void add_implicit_hydrogens();

            /**
             * @brief Set the properties of this Atom based on a .pdb format ATOM string. 
             */
            void parse_pdb(const std::string& s) override;

            /**
             * @brief Create a .pdb format string representation of this Atom. 
             */
            std::string as_pdb() const override;

            /**
             * @brief Check if this is a water molecule by checking the residue name. Both HOH and SOL are considered water molecules.
             */
            virtual bool is_water() const;

            /**
             * @brief Set the atomic element for this atom. Any spaces are removed. 
             */
            void set_element(constants::atom_t element);

            /**
             * @brief Set the atomic element for this atom. Any spaces are removed. 
             */
            void set_element(const std::string& element);

            /**
             * @brief Get the coordinates of this Atom.
             */
            Vector3<double>& coordinates();

            /**
             * @brief Get the coordinates of this Atom.
             */
            const Vector3<double>& coordinates() const;

            /**
             * @brief Get the record name of this Atom. 
             *        This is the first 6 characters of the line in the PDB file.
             *        Example: ATOM
             */
            virtual std::string get_recName() const;

            /**
             * @brief Get the mass of this atom. 
             *        If implicit hydrogens are enabled, the mass of the atom
             *        is calculated as the mass of the nucleus plus the mass of the attached hydrogen atoms.
             */
            virtual double get_mass() const;

            /**
             * @brief Get the number of protons in this atom.
             */
            unsigned int Z() const;

            /**
             * @brief Comparison function to allow this class to be a map key. 
             * 
             * @param rhs Atom to compare against.
             */
            bool operator<(const PDBAtom& rhs) const;

            /**
             * @brief Equality operator to determine if two atoms are equal.
             *        Note that this compares their unique object identifier which is generated at object creation, completely disregarding
             *        their contents. Unless a deliberate attempt at desyncing the id from the contents were made, equality of content follows
             *        from equality of id. 
             * @param rhs Atom to compare against. 
             */
            bool operator==(const PDBAtom& rhs) const;

            /**
             * @brief Equality operator to determine if two atoms are equal.
             *        Note that this is a @a content comparator, and thus determines if two atoms are equal based on their contents. 
             * @param rhs Atom to compare against. 
             */
            bool equals_content(const PDBAtom& rhs) const;

            /**
             * @brief Inequality operator to determine if two atoms are not equal.
             *        Note that this compares their unique object identifier which is generated at object creation, completely disregarding
             *        their contents. Unless a deliberate attempt at desyncing the id from the contents were made, inequality of content follows
             *        from inequality of id. 
             * @param rhs Atom to compare against. 
             */
            bool operator!=(const PDBAtom& rhs) const {return !operator==(rhs);}

            // properties as defined in https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf, page 180.
            Vector3<double> coords = {0, 0, 0};
            std::string name, altLoc, resName, iCode, charge;
            char chainID = ' ';
            constants::atom_t element = constants::atom_t::unknown;
            constants::atomic_group_t atomic_group = constants::atomic_group_t::unknown;
            double occupancy = -1, tempFactor = -1;
            int serial = -1, resSeq = -1; 

            // other properties
            double effective_charge = -1;
            int uid = -1;

        private: 
            static inline int uid_counter = 0; // global counter for unique ids
    };
}