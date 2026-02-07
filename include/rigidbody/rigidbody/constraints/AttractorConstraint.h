// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

#include <string>

namespace ausaxs::rigidbody::constraints {
    class AttractorConstraint : public DistanceConstraint {
        public: 
            /**
             * @brief Create a new constraint between a pair of atoms in the two bodies.
             * 
             * Complexity: O(n)
             */
            AttractorConstraint(observer_ptr<const data::Molecule> molecule, unsigned int ibody1, unsigned int ibody2, double target_distance);
            virtual ~AttractorConstraint() override = default;

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;

            /**
             * @brief Check if a constraint is identical to this object. 
             * 
             * @param constraint The constraint to be checked for equality. 
             */
            bool operator==(const AttractorConstraint& constraint) const;

            /**
             * @brief Generate a string representation of this constraint.
             */
            std::string to_string() const;

        private: 
            static double transform(double distance, double r_base);
    };
}