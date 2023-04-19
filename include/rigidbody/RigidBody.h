#pragma once

#include <data/Protein.h>
#include <rigidbody/constraints/Constraint.h>
#include <rigidbody/selection/BodySelectStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/constraints/ConstraintManager.h>
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
			 * @brief Apply a calibration to this rigid body. 
			 * 
			 * This will fix the solvent scattering density to the fitted value.
			 */
			void apply_calibration(std::shared_ptr<fitter::Fit> calibration);

			/**
			 * @brief Update the fitter with the current rigid body parameters.
			 */
			void update_fitter(std::shared_ptr<fitter::LinearFitter> fitter);

			std::shared_ptr<ConstraintManager> constraints;
		protected:
			std::shared_ptr<fitter::Fit> calibration = nullptr;
			std::unique_ptr<BodySelectStrategy> body_selector;
			std::unique_ptr<TransformStrategy> transform;
			std::unique_ptr<ParameterGenerationStrategy> parameter_generator;

			/**
			 * @brief Prepare the fitter for this rigidbody.
			 */
			std::shared_ptr<fitter::LinearFitter> prepare_fitter(std::string measurement_path); 

			/**
			 * @brief Small initialization function.
			 */
			void setup();
	};
}