#pragma once

#include <fitter/Fitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <utility/observer_ptr.h>

#include <vector>
#include <memory>

namespace ausaxs::fitter {
    template<typename C>
    concept fitter_t = std::is_base_of_v<Fitter, C>;

    /**
     * @brief Fit an intensity curve to a dataset. 
     * 
     * Extends a fitter with the ability to add constraints to the optimization.
     * Note that the constraint manager must manually be set with the set_constraint_manager method.
     */
    template<fitter_t T>
    class ConstrainedFitter : public T {
        public: 
            using T::T;
            ConstrainedFitter(ConstrainedFitter<T>&& other);

            [[nodiscard]] double chi2(const std::vector<double>& params) override;

            /**
             * @brief Set the constraints for the fitter. 
             * 
             * Overwrites any existing constraints.
             */
            void set_constraint_manager(std::shared_ptr<rigidbody::constraints::ConstraintManager> constraints);

            observer_ptr<rigidbody::constraints::ConstraintManager> get_constraint_manager(); 

        private: 
            std::shared_ptr<rigidbody::constraints::ConstraintManager> constraints = nullptr;
    };
}

#include <rigidbody/constraints/ConstrainedFitter.tpp>