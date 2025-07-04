// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/Constraint.h>
#include <data/DataFwd.h>

#include <iostream>

namespace ausaxs::rigidbody::constraints {
    /**
     * This class is the glue that keeps separate bodies together during the optimization. Each constraint is between two individual atoms of two different bodies, and works by 
     * adding a new term to the chi-square which is a function of the distance between the two atoms. Thus the optimizer is penalized depending on how much it separates the two 
     * constrained atoms. 
     */
    class DistanceConstraint : public Constraint {
        public: 
            /**
             * @brief Create a new constraint between a pair of atoms.
             * 
             * We use indexes since the bodies and atoms may change during the optimization.
             * Thus to use a constraint, the order of the bodies and atoms cannot change. 
             * 
             * Complexity: O(1)
             * 
             * @param protein The protein this constraint belongs to.
             * @param ibody1 The index of the first body.
             * @param ibody2 The index of the second body.
             * @param iatom1 The index of the first atom in the first body.
             * @param iatom2 The index of the second atom in the second body.
             */
            DistanceConstraint(data::Molecule* protein, unsigned int ibody1, unsigned int ibody2, unsigned int iatom1, unsigned int iatom2);

            /**
             * @brief Create a new constraint between a pair of atoms in the two bodies.
             * 
             * Complexity: O(n)
             * 
             * @param body1 The first body.
             * @param body2 The second body.
             * @param center_mass Create a constraint from the center-masses of the two bodies. If false, the two closest carbon atoms are used instead. 
             */
            DistanceConstraint(data::Molecule* protein, unsigned int ibody1, unsigned int ibody2, bool center_mass = false);

            /**
             * @brief Create a new constraint between a pair of atoms.
             * 
             * Complexity: O(n) where n is the number of atoms in the protein.
             * 
             * @param protein The protein this constraint belongs to.
             * @param atom1 The first atom.
             * @param atom2 The second atom.
             */
            DistanceConstraint(data::Molecule* protein, const data::AtomFF& atom1, const data::AtomFF& atom2);

            virtual ~DistanceConstraint() override = default;

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;

            /**
             * @brief Get the first atom of this constraint. 
             */
            const data::AtomFF& get_atom1() const;

            /**
             * @brief Get the second atom of this constraint. 
             */
            const data::AtomFF& get_atom2() const;

            /**
             * @brief Get the first body of this constraint. 
             */
            const data::Body& get_body1() const;

            /**
             * @brief Get the first body of this constraint. 
             */
            data::Body& get_body1();

            /**
             * @brief Get the second body of this constraint. 
             */
            const data::Body& get_body2() const;

            /**
             * @brief Get the second body of this constraint. 
             */
            data::Body& get_body2();

            /**
             * @brief Check if a constraint is identical to this object. 
             * 
             * @param constraint The constraint to be checked for equality. 
             */
            bool operator==(const DistanceConstraint& constraint) const;

            /**
             * @brief Generate a string representation of this constraint.
             */
            std::string to_string() const;

            friend std::ostream& operator<<(std::ostream& os, const DistanceConstraint& constraint) {os << constraint.to_string(); return os;}

            double r_base;              // The normal distance between the two atoms. 
            data::Molecule* protein;    // The protein this constraint belongs to.
            unsigned int ibody1 = -1;   // The index of the first body.
            unsigned int ibody2 = -1;   // The index of the second body.
            unsigned int iatom1 = -1;   // The index of the first atom.
            unsigned int iatom2 = -1;   // The index of the second atom.

        private: 
            struct AtomLoc {int body, atom;};

            /**
             * @brief Transforms a distance into a proper constraint for least-squares fitting. 
             * 
             * @param offset The radial offset between the new and original positions. 
             */
            static double transform(double offset);

            /**
             * @brief Find the bodies containing the argument atoms.
             *        This is linear in the total number of atoms. 
             */
            std::pair<AtomLoc, AtomLoc> find_host_bodies(const data::AtomFF& atom1, const data::AtomFF& atom2) const;
    };
}