#pragma once

#include <rigidbody/constraints/Constraint.h>
#include <hist/Histogram.h>

class Protein;
namespace rigidbody {

    /**
     * @brief Overlap constraint. 
     * 
     * This constraint will try to reduce the overlap between atoms in different bodies. 
     * More specifically an exponentially decaying function is used as a weight. The product of this weight and the initial distance between the atoms is then used as a target. 
     * The squared deviation from this target is the chi2 contribution of this constraint.
     */
    class OverlapConstraint : public Constraint {
        public:
            OverlapConstraint() = default;

            OverlapConstraint(Protein* protein);

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;

        private: 
            Protein* protein;
            hist::Histogram target;
            hist::Histogram weights;

            /**
             * @brief Initialize the target distribution.
             */
            void initialize();
    };
}