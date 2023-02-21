#pragma once

#include <fitter/HydrationFitter.h>
#include <rigidbody/Constraint.h>

namespace fitter {
    /**
     * @brief Fit an intensity curve to a dataset. 
     * 
     * Three parameters will be fitted: 
     *    a: The slope of the curve.
     *    b: The intercept of the curve.
     *    c: The scattering length of the hydration shell.
     * 
     * This fitter also supports atomic constraints. 
     */
    class ConstrainedFitter : public HydrationFitter {
        public: 
            using HydrationFitter::HydrationFitter;

            [[nodiscard]] double chi2(const std::vector<double>& params) override;

            /**
             * @brief Add a constraint to the fitter. 
             */
            void add_constraint(const rigidbody::Constraint& constraint);

            /**
             * @brief Add a constraint to the fitter. 
             */
            void add_constraint(rigidbody::Constraint&& constraint);

            /**
             * @brief Set the constraints for the fitter. 
             * 
             * Overwrites any existing constraints.
             */
            void set_constraints(std::vector<rigidbody::Constraint>&& constraints);

            /**
             * @brief Set the constraints for the fitter. 
             * 
             * Overwrites any existing constraints.
             */
            void set_constraints(const std::vector<rigidbody::Constraint>& constraints);

            /**
             * @brief Get the constraints. 
             * 
             * @return The constraints. 
             */
            [[nodiscard]] const std::vector<rigidbody::Constraint>& get_constraints() const;            

        private: 
            std::vector<rigidbody::Constraint> constraints;
    };
}
