#pragma once

#include <memory>

#include <data/Protein.h>
#include <rigidbody/Constraint.h>
#include <rigidbody/BodySelectStrategy.h>
#include <rigidbody/TransformationStrategy.h>
#include <rigidbody/ParameterGenerationStrategy.h>
#include <fitter/HydrationFitter.h>

namespace rigidbody {
	class RigidBody {
		public:
			/**
			 * @brief Construtor. 
			 * 
			 * Prepare a new rigid body for optimization. 
			 * 
			 * @param protein The protein to be optimized. 
			 */
			explicit RigidBody(Protein& protein);

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

			Protein& protein;
			std::vector<Constraint> constraints;

		private:
			std::unique_ptr<BodySelectStrategy> body_selector;
			std::unique_ptr<TransformationStrategy> transform;
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