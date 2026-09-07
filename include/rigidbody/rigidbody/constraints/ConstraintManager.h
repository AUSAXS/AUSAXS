// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>
#include <rigidbody/constraints/OverlapConstraint.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/RigidbodyFwd.h>
#include <settings/RigidBodySettings.h>
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
         * @brief Generate automatic constraints using the given strategy.
         */
        void generate_constraints(settings::rigidbody::ConstraintGenerationStrategyChoice choice);

        /**
         * @brief Generate automatic constraints using a custom generator, appending them to whatever constraints are already present.
         *
         * Generated constraints are appended rather than assigned over the list, so a constraint the user declared beforehand survives - deleting it would
         * throw away a real chi2 tether and a real transform-propagation link, and the declaration and the generation are equally deliberate statements.
         * Appending is only correct because generation happens at most once; a second run is rejected rather than silently doubling every bond it finds.
         *
         * @throws If constraints have already been generated.
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

        /**
         * @brief Mark the per-body map as stale, so the next access rebuilds it.
         *
         * Must be called by anything that adds or removes a body. The map is keyed by body index, so a molecule that gained a body would otherwise be missing
         * an entry for it, and one that lost a body would keep a stale entry. Pose changes need no invalidation: the map records constraint membership, not
         * where the bodies are.
         */
        void invalidate() {dirty = true;}

        observer_ptr<const data::Molecule> molecule;
        std::vector<std::unique_ptr<IDistanceConstraint>> discoverable_constraints;
        std::vector<std::unique_ptr<Constraint>> non_discoverable_constraints;

        private:
            /**
            * @brief Rebuild the per-body constraint map, if it has been invalidated since the last access.
            *
            * This map allows us to quickly find all constraints that apply to a given body without having to iterate over all constraints.
            */
            void refresh() const;

            // Whether a generator has already run. Generation is additive, so allowing it twice would duplicate every constraint the first run found.
            bool generated = false;

            // List of all discoverable constraints associated with each body. It is a cache of what the constraint list and molecule already say, so reading it
            // is a const operation even when it has to be rebuilt first.
            mutable std::unordered_map<unsigned int, std::vector<observer_ptr<IDistanceConstraint>>> distance_constraints_map;
            mutable bool dirty = true;
    };
}