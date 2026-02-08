// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <fitter/Fitter.h>
#include <fitter/SmartFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <utility/observer_ptr.h>

#include <cassert>
#include <vector>

namespace ausaxs::fitter {
    /**
     * @brief Fit an intensity curve to a dataset. 
     * 
     * Extends a fitter with the ability to add constraints to the optimization.
     * Note that the constraint manager must manually be set with the set_constraint_manager method.
     */
    struct ConstrainedFitter : SmartFitter {
        // warn of pointer issues with the embedded function
        ConstrainedFitter(ConstrainedFitter&&) noexcept = delete;
        ConstrainedFitter(const ConstrainedFitter&) = delete;
        ConstrainedFitter& operator=(ConstrainedFitter&& other) noexcept = delete;
        ConstrainedFitter& operator=(const ConstrainedFitter& other) = delete;

        template <typename... Args, typename = std::enable_if_t<std::is_constructible_v<SmartFitter, Args...>>>
        ConstrainedFitter(observer_ptr<rigidbody::constraints::ConstraintManager> constraints, Args&&... args) : SmartFitter(std::forward<Args>(args)...), constraints(constraints) {
            assert(constraints != nullptr && "ConstrainedFitter: Constraint manager must not be null.");
        }

        [[nodiscard]] double chi2(const std::vector<double>& params) override {
            return SmartFitter::chi2(params) + constraints_chi2();
        }

        /**
         * @brief Perform a chi2-only fit with cached constraint evaluation.
         * 
         * Since constraints are geometric (independent of fitted parameters a, b, cw),
         * we evaluate them once and cache the result to avoid redundant computation
         * during the parameter optimization.
         */
        [[nodiscard]] double fit_chi2_only() override {
            double cached_constraint_value = constraints->evaluate();
            constraints_chi2 = [cached_constraint_value] () {return cached_constraint_value;}; // use cached value
            double result = SmartFitter::fit_chi2_only();
            constraints_chi2 = [this] () {return constraints->evaluate();}; // re-enable dynamic evaluation
            return result;
        }

        std::unique_ptr<FitResult> unconstrained_fit() {
            constraints_chi2 = [] () {return 0;}; // disable constraints for the unconstrained fit
            auto res = SmartFitter::fit();
            constraints_chi2 = [this] () {return constraints->evaluate();}; // re-enable constraints
            return res;
        }

        observer_ptr<rigidbody::constraints::ConstraintManager> constraints;
        std::function<double()> constraints_chi2 = [this] () {return constraints->evaluate();};
    };
}