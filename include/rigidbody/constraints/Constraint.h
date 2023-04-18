#pragma once

namespace rigidbody {
    /**
     * @brief Constraint. 
     * 
     * A constraint is a function that is evaluated by the optimizer in each step. The optimizer will try to minimize the chi2 of the system, and thus the constraints will
     * be minimized as well.
     */
    class Constraint {
        public: 
            /**
             * @brief Default constructor.
             */
            Constraint() {}

            /**
             * @brief Evaluate this constraint. This method is called by the optimizer in each step to evaluate the chi2 contribution of this constraint.
             */
            virtual double evaluate() const = 0;
    };
}