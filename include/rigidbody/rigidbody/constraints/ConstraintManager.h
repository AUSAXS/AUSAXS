// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <data/DataFwd.h>

#include <vector>
#include <unordered_map>

namespace ausaxs::rigidbody::constraints {
    struct ConstraintManager {
        ConstraintManager(observer_ptr<const data::Molecule> molecule);
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
         * @brief Add a new overlap constraint.
         */
        void add_constraint(OverlapConstraint&& constraint);

        /**
         * @brief Add a new overlap constraint.
         */
        void add_constraint(const OverlapConstraint& constraint);

        /**
         * @brief Add a new distance constraint.
         */
        void add_constraint(DistanceConstraint&& constraint);

        /**
         * @brief Add a new distance constraint.
         */
        void add_constraint(const DistanceConstraint& constraint);

        /**
         * @brief Evaluate all constraints.
         * 
         * @return The chi2 contribution of all constraints.
         */
        double evaluate() const;

        observer_ptr<const data::Molecule> molecule;
        OverlapConstraint overlap_constraint;
        std::vector<DistanceConstraint> distance_constraints;
        std::unordered_map<unsigned int, std::vector<std::reference_wrapper<DistanceConstraint>>> distance_constraints_map; // Maps a body index to all its constraints

        private:
            /**
            * @brief Generate a map of constraints for each body.
            * 
            * This map allows us to quickly find all constraints that apply to a given body without having to iterate over all constraints.
            */
            void update_constraint_map();
    };
}