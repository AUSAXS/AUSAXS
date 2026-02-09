// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/IDistanceConstraint.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

#include <utility>

namespace ausaxs::rigidbody::constraints {
    /**
     * @brief Distance constraint between individual atoms with symmetry support.
     * 
     * This constraint maintains a target distance between specific atoms of different bodies,
     * accounting for symmetry transformations that create virtual bonds between symmetric partners.
     */
    class DistanceConstraintAtom : public IDistanceConstraint {
        public:
            /**
             * @brief Create an atom distance constraint between a pair of atoms with symmetry support.
             * 
             * @param molecule The molecule this constraint belongs to.
             * @param ibody1 Index of the first body.
             * @param ibody2 Index of the second body.
             * @param iatom1 Index of the first atom in body1.
             * @param iatom2 Index of the second atom in body2.
             * @param isym1 Symmetry index for body1 (default: no symmetry).
             * @param isym2 Symmetry index for body2 (default: no symmetry).
             */
            DistanceConstraintAtom(
                observer_ptr<const data::Molecule> molecule, 
                int ibody1, int iatom1, int ibody2, int iatom2,
                std::pair<int, int> isym1 = {-1, -1}, std::pair<int, int> isym2 = {-1, -1}
            );

            DistanceConstraintAtom(
                observer_ptr<const data::Molecule> molecule, 
                int ibody1, int ibody2,
                std::pair<int, int> isym1 = {-1, -1}, std::pair<int, int> isym2 = {-1, -1}
            );

            /**
             * @brief Create a constraint between two specific atoms by reference.
             * 
             * This constructor searches through all bodies to find which ones contain the given atoms.
             * 
             * @param molecule The molecule this constraint belongs to.
             * @param atom1 The first atom.
             * @param atom2 The second atom.
             */
            DistanceConstraintAtom(
                observer_ptr<const data::Molecule> molecule,
                const data::AtomFF& atom1,
                const data::AtomFF& atom2
            );

            virtual ~DistanceConstraintAtom() override = default;

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;

        protected:
            /**
             * @brief Calculate the current distance between the (symmetry-transformed) atoms.
             */
            double evaluate_distance(int iatom1, int iatom2) const;
    };
}