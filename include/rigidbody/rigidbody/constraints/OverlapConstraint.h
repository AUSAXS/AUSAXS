#pragma once

#include <rigidbody/constraints/Constraint.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <utility/observer_ptr.h>
#include <data/DataFwd.h>

namespace rigidbody::constraints {

    /**
     * @brief Overlap constraint. 
     * 
     * This constraint will try to reduce the overlap between atoms in different bodies. 
     * More specifically an exponentially decaying function is used as a weight. The product of this weight and the initial distance between the atoms is then used as a target. 
     * The squared deviation from this target is the chi2 contribution of this constraint.
     */
    class OverlapConstraint : public Constraint {
        public:
            OverlapConstraint(data::Molecule* protein);

            virtual ~OverlapConstraint() override;

            /**
             * @brief Evaluate this constraint for the current positions. 
             * 
             * @return The chi2 contribution of this constraint.
             */
            double evaluate() const override;

            bool operator==(const OverlapConstraint& other) const;

            static void set_overlap_function(std::function<double(double)> func);

        protected:
            static double weight(double r);

        private: 
            inline static std::function<double(double)> overlap_function = [](double r) {return std::exp(-5*r);};
            observer_ptr<data::Molecule> protein;
            std::vector<double> target;
            std::vector<double> weights;
            std::vector<double> axis;

            /**
             * @brief Initialize the target distribution.
             */
            void initialize();
    };
}