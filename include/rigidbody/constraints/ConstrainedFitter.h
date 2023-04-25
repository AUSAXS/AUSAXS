#pragma once

#include <fitter/HydrationFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>

namespace rigidbody {
    class DistanceConstraint;
}

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
             * @brief Set the constraints for the fitter. 
             * 
             * Overwrites any existing constraints.
             */
            void set_constraint_manager(std::shared_ptr<rigidbody::ConstraintManager> constraints);

        private: 
            std::shared_ptr<rigidbody::ConstraintManager> constraints;
    };
}

#include <rigidbody/constraints/ConstrainedFitter.tpp>