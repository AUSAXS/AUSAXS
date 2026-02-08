// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::constraints {
    class AttractorConstraint : public DistanceConstraintCM {
        public: 
            /**
             * @brief Create a new constraint between a pair of atoms in the two bodies.
             * 
             * Complexity: O(n)
             */
            AttractorConstraint(
                observer_ptr<const data::Molecule> molecule, double target_distance, 
                int ibody1, int ibody2, std::pair<int, int> isym1 = {-1, -1}, std::pair<int, int> isym2 = {-1, -1}
            );
            virtual ~AttractorConstraint() override = default;

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;

        private: 
            static double transform(double distance, double r_base);
    };
}