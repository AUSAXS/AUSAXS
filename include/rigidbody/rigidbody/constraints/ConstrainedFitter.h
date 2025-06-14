// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <fitter/Fitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <utility/observer_ptr.h>

#include <vector>

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
    struct ConstrainedFitter : T {
        template <typename... Args, typename = std::enable_if_t<std::is_constructible_v<T, Args...>>>
        ConstrainedFitter(observer_ptr<rigidbody::constraints::ConstraintManager> constraints, Args&&... args) : T(std::forward<Args>(args)...), constraints(constraints) {
            assert(constraints != nullptr && "ConstrainedFitter: Constraint manager must not be null.");
        }

        ConstrainedFitter(ConstrainedFitter<T>&& other) : T(std::move(other)), constraints(other.constraints) {}

        [[nodiscard]] double chi2(const std::vector<double>& params) override {
            return T::chi2(params) + constraints->evaluate();
        }
        observer_ptr<rigidbody::constraints::ConstraintManager> constraints;
    };
}