#pragma once

#include <fitter/HydrationFitter.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>

namespace fitter {
    template<typename C>
    concept fitter_t = std::is_base_of_v<Fitter, C>;

    /**
     * @brief Fit an intensity curve to a dataset. 
     * 
     * Extends a fitter with the ability to add constraints to the optimization.
     */
    template<fitter_t T>
    class ConstrainedFitter : public T {
        public: 
            using T::T;
            [[nodiscard]] double chi2(const std::vector<double>& params) override;

            /**
             * @brief Add a constraint to the fitter. 
             */
            void add_constraint(std::shared_ptr<rigidbody::DistanceConstraint> constraint);

            /**
             * @brief Add a constraint to the fitter. 
             */
            void add_constraint(std::shared_ptr<rigidbody::DistanceConstraint> constraint);

            /**
             * @brief Add a constraint to the fitter. 
             */
            void add_constraint(rigidbody::DistanceConstraint&& constraint);

            /**
             * @brief Set the constraints for the fitter. 
             * 
             * Overwrites any existing constraints.
             */
            void set_constraints(std::vector<std::shared_ptr<rigidbody::DistanceConstraint>> constraints);

            /**
             * @brief Set the constraints for the fitter. 
             * 
             * Overwrites any existing constraints.
             */
            void set_constraints(std::vector<std::shared_ptr<rigidbody::DistanceDistanceConstraint>> constraints);

            /**
             * @brief Get the constraints. 
             * 
             * @return The constraints. 
             */
            [[nodiscard]] const std::vector<std::shared_ptr<rigidbody::DistanceConstraint>>& get_constraints() const;

        private: 
            std::vector<std::shared_ptr<rigidbody::DistanceConstraint>> distance_constraints;
            std::unique_ptr<rigidbody::OverlapConstraint> overlap_constraint;
    };
}

#include <rigidbody/constraints/ConstrainedFitter.tpp>