// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/RigidbodyFwd.h>
#include <data/DataFwd.h>

#include <memory>
#include <vector>
#include <unordered_map>

namespace ausaxs::rigidbody::constraints {
    /**
     * @brief Owns and evaluates the constraints applied during a rigidbody optimization.
     *
     * Constraints are kept in two groups. Discoverable constraints are distance constraints tied
     * to specific bodies; they are additionally indexed in a per-body map so that all constraints
     * affecting a given body can be looked up cheaply without scanning the full list. Non-
     * discoverable constraints (such as the overlap constraint) are not tied to a body pair and
     * are only evaluated globally. evaluate() sums the chi2 penalty contributed by both groups.
     */
    struct ConstraintManager {
        ConstraintManager(observer_ptr<const Rigidbody> rigidbody);
        ~ConstraintManager();

        /**
         * @brief Generate automatic constraints based on the currently selected constraint generation strategy in the settings.
         */
        void generate_constraints();

        /**
         * @brief Generate automatic constraints using a custom generator.
         */
        void generate_constraints(std::unique_ptr<ConstraintGenerationStrategy> generator);

        /**
         * @brief Add a new constraint.
         */
        void add_constraint(std::unique_ptr<Constraint> constraint);

        /**
         * @brief Evaluate all constraints.
         * 
         * @return The chi2 contribution of all constraints.
         */
        double evaluate() const;

        /**
         * @brief Get all discoverable distance constraints that involve the given body.
         */
        const std::vector<observer_ptr<IDistanceConstraint>>& get_body_constraints(unsigned int ibody) const;

        observer_ptr<const data::Molecule> molecule;
        std::vector<std::unique_ptr<IDistanceConstraint>> discoverable_constraints;
        std::vector<std::unique_ptr<Constraint>> non_discoverable_constraints;

        private:
            // List of all discoverable constraints associated with each body. 
            std::unordered_map<unsigned int, std::vector<observer_ptr<IDistanceConstraint>>> distance_constraints_map;

            /**
            * @brief Generate a map of constraints for each body.
            * 
            * This map allows us to quickly find all constraints that apply to a given body without having to iterate over all constraints.
            */
            void update_constraint_map();
    };
}