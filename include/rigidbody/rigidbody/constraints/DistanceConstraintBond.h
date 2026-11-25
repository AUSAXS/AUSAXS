// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::constraints {
    class DistanceConstraintBond : public DistanceConstraintAtom {
        public: 
            /**
             * @brief Create a backbone bond constraint between the two bodies.
             *
             * The representative C-alpha pair is selected automatically: only the terminal C-alphas of each body are eligible, and of those
             * only pairs sharing a chain and lying within a few residues of each other in sequence. Residues missing from the model are
             * allowed for, so two fragments separated by a disordered stretch still bond.
             * Throws if no suitable pair is found or if the selected atoms are too far apart for the number of residues between them.
             */
            DistanceConstraintBond(observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2);

            ~DistanceConstraintBond() override = default;

            /**
             * @brief Check whether a backbone bond can be formed between the two bodies.
             *
             * Equivalent to the constructor succeeding, but without throwing, so many pairs can be probed cheaply.
             */
            static bool can_bond(observer_ptr<const data::Molecule> molecule, int ibody1, int ibody2);

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;
    };
}