#pragma once

#include <fitter/Fitter.h>
#include <rigidbody/constraints/ConstraintManager.h>

#include <vector>
#include <memory>
#include <concepts>

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

            rigidbody::ConstraintManager* get_constraint_manager(); 

        private: 
            std::shared_ptr<rigidbody::ConstraintManager> constraints = std::make_unique<rigidbody::ConstraintManager>();
    };
}

#include <rigidbody/constraints/ConstrainedFitter.tpp>