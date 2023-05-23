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

            // /**
            //  * @brief Add a constraint to this manager.
            //  */
            // template<typename T, typename = std::enable_if_t<std::is_same_v<T, DistanceConstraint>>>
            // void add_constraint(T&& constraint) {
            //     distance_constraints.push_back(std::forward(constraint));
            //     generate_constraint_map();
            // }

            template<typename T, typename = std::enable_if_t<std::is_base_of_v<T, Constraint>>>
            void add_constraint(T&& constraint);

            // /**
            //  * @brief Add a constraint to this manager.
            //  */
            // template<typename T, typename = std::enable_if_t<std::is_same_v<T, OverlapConstraint>>>
            // void add_constraint(T&& constraint) {
            //     overlap_constraint = std::forward(constraint);
            // }

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