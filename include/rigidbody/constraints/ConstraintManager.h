#pragma once

#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/OverlapConstraint.h>

#include <vector>
#include <memory>
#include <unordered_map>

class Protein;
namespace rigidbody {
    class ConstraintManager {
        public:
            ConstraintManager();

            /**
             * @brief Construct a new Constraint Manager for a given protein.
             */
            ConstraintManager(Protein* protein);

            ~ConstraintManager();

            void add_constraint(OverlapConstraint&& constraint);
            void add_constraint(const OverlapConstraint& constraint);
            void add_constraint(DistanceConstraint&& constraint);
            void add_constraint(const DistanceConstraint& constraint);

            /**
             * @brief Evaluate all constraints.
             * 
             * @return The chi2 contribution of all constraints.
             */
            double evaluate() const;

            Protein* protein;
            OverlapConstraint overlap_constraint;                                                        // The overlap constraint
            std::vector<DistanceConstraint> distance_constraints;                                        // All distance constraints
			std::unordered_map<unsigned int, std::vector<DistanceConstraint*>> distance_constraints_map; // Maps a body index to all its constraints

        private:
            /**
			 * @brief Generate a map of constraints for each body.
			 * 
			 * This map allows us to quickly find all constraints that apply to a given body without having to iterate over all constraints.
			 */
            void generate_constraint_map();
    };
}