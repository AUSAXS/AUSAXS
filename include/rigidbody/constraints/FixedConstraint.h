#pragma once

#include <rigidbody/constraints/Constraint.h>
#include <utility/observer_ptr.h>
#include <data/DataFwd.h>

namespace rigidbody::constraints {
    /**
     * This constraint can be used to link different bodies such that any transformation applied to one will also be applied to the other.
     * This is principially equivalent to the DistanceConstraint, except this constraint does not affect the chi2. 
     */
    class FixedConstraint : public Constraint {
        public: 
            /**
             * @brief Create a new constraint between the center of mass of two bodies.
             * 
             * Complexity: O(1)
             * 
             * @param protein The protein this constraint belongs to.
             * @param ibody1 The index of the first body.
             * @param ibody2 The index of the second body.
             */
            FixedConstraint(observer_ptr<data::Molecule> protein, unsigned int ibody1, unsigned int ibody2);

            virtual ~FixedConstraint() override = default;

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;

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

            observer_ptr<data::Molecule> protein;   // The protein this constraint belongs to.
            unsigned int ibody1, ibody2;            // The indices of the bodies.
    };
}