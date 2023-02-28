#pragma once

#include <data/Protein.h>
#include <rigidbody/Constraint.h>
#include <rigidbody/selection/BodySelectStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <fitter/HydrationFitter.h>

#include <memory>
#include <unordered_map>

namespace rigidbody {
	class RigidBody : public Protein {
		public:
			RigidBody(Protein&& protein);

			RigidBody(const Protein& protein);

			/**
			 * @brief Perform a rigid-body optimization for this structure. 
			 */
			std::shared_ptr<fitter::Fit> optimize(std::string measurement_path);

			/**
			 * @brief Add a constraint to this rigid body. 
			 */
			void add_constraint(std::shared_ptr<rigidbody::Constraint> constraint);
			
			/**
			 * @brief Add a constraint to this rigid body. 
			 */
			void add_constraint(rigidbody::Constraint&& constraint);

			/**
			 * @brief Add a constraint to this rigid body. 
			 */
			void add_constraint(unsigned int ibody1, unsigned int ibody2, unsigned int iatom1, unsigned int iatom2);

			/**
			 * @brief Generate a set of simple constraints. 
			 * 
			 * This method will generate and add one constraint between each pair of connected bodies. 
			 */
			void generate_simple_constraints();

			/**
			 * @brief Generate a map of constraints for each body.
			 * 
			 * This map allows us to quickly find all constraints that apply to a given body without having to iterate over all constraints.
			 */
            void generate_constraint_map();

			/**
			 * @brief Apply a calibration to this rigid body. 
			 * 
			 * This will fix the solvent scattering density to the fitted value.
			 */
			void apply_calibration(std::shared_ptr<fitter::Fit> calibration);

			/**
			 * @brief Update the fitter with the current rigid body parameters.
			 */
			void update_fitter(std::shared_ptr<fitter::LinearFitter> fitter);

			std::vector<std::shared_ptr<Constraint>> get_constraints() const;
			std::shared_ptr<Constraint> get_constraint(unsigned int index) const;

			std::unordered_map<unsigned int, std::vector<std::shared_ptr<Constraint>>> constraint_map;
		protected:
			std::shared_ptr<fitter::Fit> calibration = nullptr;
			std::unique_ptr<BodySelectStrategy> body_selector;
			std::unique_ptr<TransformStrategy> transform;
			std::unique_ptr<ParameterGenerationStrategy> parameter_generator;
			std::vector<std::shared_ptr<Constraint>> constraints;

			/**
			 * @brief Prepare the fitter for this rigidbody.
			 */
			std::shared_ptr<fitter::LinearFitter> prepare_fitter(std::string measurement_path); 

			/**
			 * @brief Perform a single step of the optimization, and calculate the resulting chi2 value. 
			 */
			double chi2(fitter::HydrationFitter& fitter) const;

			/**
			 * @brief Small initialization function.
			 */
			void setup();
	};
}