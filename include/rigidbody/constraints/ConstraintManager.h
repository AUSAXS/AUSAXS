#pragma once

#include <vector>
#include <memory>
#include <unordered_map>

class Protein;
namespace rigidbody {
    class OverlapConstraint;
    class DistanceConstraint;
    class ConstraintManager {
        public:
            /**
             * @brief Construct a new Constraint Manager for a given protein.
             */
            ConstraintManager(Protein* protein);

            /**
             * @brief Add a constraint to this manager.
             */
            void add_constraint(std::shared_ptr<DistanceConstraint> constraint);

            /**
             * @brief Add a constraint to this manager.
             */
            void add_constraint(std::shared_ptr<OverlapConstraint> constraint);

            /**
             * @brief Evaluate all constraints.
             * 
             * @return The chi2 contribution of all constraints.
             */
            double evaluate() const;

            Protein* protein;
            std::shared_ptr<OverlapConstraint> overlap_constraint;                  // The overlap constraint
            std::vector<std::shared_ptr<DistanceConstraint>> distance_constraints;  // All distance constraints
			std::unordered_map<unsigned int, std::vector<std::shared_ptr<DistanceConstraint>>> distance_constraints_map; // Maps a body index to all its constraints

        private:
            /**
			 * @brief Generate a map of constraints for each body.
			 * 
			 * This map allows us to quickly find all constraints that apply to a given body without having to iterate over all constraints.
			 */
            void generate_constraint_map();
    };
}