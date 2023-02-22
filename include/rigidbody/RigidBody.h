#pragma once

#include <memory>

#include <data/Protein.h>
#include <rigidbody/Constraint.h>
#include <rigidbody/selection/BodySelectStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <fitter/HydrationFitter.h>

namespace rigidbody {
	class RigidBody : private Protein {
		public:
			friend TransformStrategy;
			friend BodySelectStrategy;
			friend ParameterGenerationStrategy;

			explicit RigidBody(std::string input);

			explicit RigidBody(std::vector<Body>&& bodies);

			/**
			 * @brief Perform a rigid-body optimization for this structure. 
			 */
			void optimize(std::string measurement_path);

			/**
			 * @brief Add a constraint to this rigid body. 
			 */
			void add_constraint(const rigidbody::Constraint& constraint);

			/**
			 * @brief Add a constraint to this rigid body. 
			 */
			void add_constraint(rigidbody::Constraint&& constraint);

			/**
			 * @brief Generate a set of simple constraints. 
			 * 
			 * This method will generate and add one constraint between each pair of connected bodies. 
			 */
			void generate_simple_constraints();

			std::vector<Constraint> constraints;

		private:
			std::unique_ptr<BodySelectStrategy> body_selector;
			std::unique_ptr<TransformStrategy> transform;
			std::unique_ptr<ParameterGenerationStrategy> parameter_generator;

			/**
			 * @brief Perform a single step of the optimization, and calculate the resulting chi2 value. 
			 */
			double chi2(fitter::HydrationFitter& fitter) const;

			/**
			 * @brief Rotate a body with the currently chosen transformation strategy. 
			 */
			void rotate();

			/**
			 * @brief Translate a body with the currently chosen transformation strategy. 
			 */
			void translate();
	};
}